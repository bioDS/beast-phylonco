package phylonco.beast.evolution.readcountmodel;

import java.io.IOException;
import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import mutablealignment.MutableAlignment;

/**
 * MCMC that incrementally adds alignment sites during the run.
 * Starts with a subset of sites (patternWeight=1), remaining sites have patternWeight=0.
 * At scheduled intervals, a batch of new sites is activated and a Gibbs site sampler
 * is called to initialise genotypes at all sites conditional on the current tree.
 *
 * <p>Example XML usage:</p>
 * <pre>
 * &lt;run id="mcmc" spec="phylonco.beast.evolution.readcountmodel.DataTemperedMCMC"
 *      chainLength="2000000"
 *      mutableAlignment="@A"
 *      gibbsSiteOperator="@GibbsSiteOp"
 *      initialSites="100"
 *      sitesPerStep="100"
 *      iterationsPerStep="50000"&gt;
 *
 *     &lt;operator id="GibbsSiteOp"
 *              spec="phylonco.beast.evolution.readcountmodel.GibbsSiteOperator"
 *              weight="1.0"
 *              sampleAllSites="true"
 *              mutableAlignment="@A"
 *              tree="@psi"
 *              siteModel="@siteModel"
 *              readCountModel="@readCountLikelihood"
 *              readCount="@readCounts"/&gt;
 *
 *     &lt;!-- other operators, state, posterior, loggers as usual --&gt;
 * &lt;/run&gt;
 * </pre>
 */
@Description("MCMC that incrementally adds alignment sites during the run. " +
        "Starts with a subset of sites (patternWeight=1), remaining sites have patternWeight=0. " +
        "At scheduled intervals, a batch of new sites is activated and a Gibbs site sampler " +
        "is called to initialise genotypes at all sites conditional on the current tree.")
public class DataTemperedMCMC extends MCMC {

    public Input<MutableAlignment> alignmentInput = new Input<>(
            "mutableAlignment",
            "the mutable alignment whose site weights will be controlled",
            Input.Validate.REQUIRED);

    public Input<GibbsSiteOperator> gibbsOperatorInput = new Input<>(
            "gibbsSiteOperator",
            "Gibbs site operator used to initialise genotypes when sites are added",
            Input.Validate.REQUIRED);

    public Input<Integer> initialSitesInput = new Input<>(
            "initialSites",
            "number of sites to activate at the start",
            100);

    public Input<Integer> sitesPerStepInput = new Input<>(
            "sitesPerStep",
            "number of sites to add at each step",
            100);

    public Input<Long> iterationsPerStepInput = new Input<>(
            "iterationsPerStep",
            "number of MCMC iterations between adding batches of sites",
            50000L);

    private MutableAlignment alignment;
    private GibbsSiteOperator gibbsOperator;
    private int totalSites;
    private int activeSites;
    private int[] siteOrder; // random permutation of site indices

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        alignment = alignmentInput.get();
        gibbsOperator = gibbsOperatorInput.get();
        totalSites = alignment.getSiteCount();

        // Create random permutation of site indices
        siteOrder = new int[totalSites];
        for (int i = 0; i < totalSites; i++) {
            siteOrder[i] = i;
        }
        shuffleArray(siteOrder);

        // Activate initial sites
        activeSites = Math.min(initialSitesInput.get(), totalSites);
        updatePatternWeights();

        Log.info.println("DataTemperedMCMC: starting with " + activeSites +
                " of " + totalSites + " sites active");
    }

    @Override
    protected void doLoop() throws IOException {
        int corrections = 0;
        final boolean isStochastic = posterior.isStochastic();

        if (burnIn > 0) {
            Log.warning.println("Please wait while BEAST takes " + burnIn + " pre-burnin samples");
        }

        long nextStepIteration = iterationsPerStepInput.get();

        for (long sampleNr = -burnIn; sampleNr <= chainLength; sampleNr++) {

            // Check if it's time to add more sites
            if (sampleNr > 0 && sampleNr >= nextStepIteration && activeSites < totalSites) {
                int oldActive = activeSites;
                activeSites = Math.min(activeSites + sitesPerStepInput.get(), totalSites);
                updatePatternWeights();

                Log.info.println("DataTemperedMCMC: iteration " + sampleNr +
                        " — activated sites " + oldActive + " -> " + activeSites +
                        " / " + totalSites);

                // Update cached cumulative weights then Gibbs sample all active sites
                gibbsOperator.updateCumulativeWeights();
                gibbsOperator.proposal();

                // Recalculate posterior with new sites
                oldLogLikelihood = state.robustlyCalcPosterior(posterior);

                nextStepIteration = sampleNr + iterationsPerStepInput.get();
            }

            final Operator operator = propagateState(sampleNr);

            if (debugFlag && sampleNr % 1 == 0 || sampleNr % 10000 == 0) {
                final double originalLogP = isStochastic ? posterior.getNonStochasticLogP() : oldLogLikelihood;
                final double logLikelihood = isStochastic ? state.robustlyCalcNonStochasticPosterior(posterior) : state.robustlyCalcPosterior(posterior);
                if (isTooDifferent(logLikelihood, originalLogP)) {
                    reportLogLikelihoods(posterior, "");
                    Log.err.println("At sample " + sampleNr + "\nLikelihood incorrectly calculated: " + originalLogP + " != " + logLikelihood
                            + "(" + (originalLogP - logLikelihood) + ")"
                            + " Operator: " + operator.getName());
                }
                if (sampleNr > NR_OF_DEBUG_SAMPLES * 3) {
                    debugFlag = false;
                    if (isTooDifferent(logLikelihood, originalLogP)) {
                        corrections++;
                        if (corrections > 100) {
                            Log.err.println("Too many corrections. There is something seriously wrong that cannot be corrected");
                            state.storeToFile(sampleNr);
                            operatorSchedule.storeToFile();
                            System.exit(1);
                        }
                        oldLogLikelihood = state.robustlyCalcPosterior(posterior);
                    }
                } else {
                    if (isTooDifferent(logLikelihood, originalLogP)) {
                        state.storeToFile(sampleNr);
                        operatorSchedule.storeToFile();
                        System.exit(1);
                    }
                }
            } else {
                if (sampleNr >= 0) {
                    operator.optimize(logAlpha);
                }
            }
            callUserFunction(sampleNr);

            if (storeEvery > 0 && (sampleNr + 1) % storeEvery == 0 || sampleNr == chainLength) {
                state.robustlyCalcNonStochasticPosterior(posterior);
                state.storeToFile(sampleNr);
                operatorSchedule.storeToFile();
            }

            if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
                throw new RuntimeException("Encountered a positive infinite posterior. This is a sign there may be numeric instability in the model.");
            }
        }
        if (corrections > 0) {
            Log.err.println("\n\nNB: " + corrections + " posterior calculation corrections were required. This analysis may not be valid!\n\n");
        }
    }

    /**
     * Set patternWeight to 1 for active sites, 0 for inactive sites.
     */
    private void updatePatternWeights() {
        int[] weights = alignment.getWeights();
        Arrays.fill(weights, 0);
        for (int i = 0; i < activeSites; i++) {
            weights[siteOrder[i]] = 1;
        }
    }

    private boolean isTooDifferent(double logLikelihood, double originalLogP) {
        return Math.abs(logLikelihood - originalLogP) > 1e-6;
    }

    /**
     * Fisher-Yates shuffle.
     */
    private static void shuffleArray(int[] array) {
        for (int i = array.length - 1; i > 0; i--) {
            int j = Randomizer.nextInt(i + 1);
            int temp = array[i];
            array[i] = array[j];
            array[j] = temp;
        }
    }
}
