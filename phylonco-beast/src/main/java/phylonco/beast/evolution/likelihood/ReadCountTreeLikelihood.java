package phylonco.beast.evolution.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;

import java.util.List;

/**
 * Tree likelihood that integrates the latent single-cell genotype alignment OUT of the
 * posterior analytically, instead of sampling it as a {@code MutableAlignment} StateNode.
 *
 * <p>Each tip's Felsenstein partial-likelihood vector is set to
 * {@code P(read counts at cell t, site j | genotype g)} for every genotype {@code g}, taken
 * from {@link LikelihoodReadCountModel#logLiklihoodRC}. Standard pruning over the GT16/GT10
 * substitution model then marginalises the genotypes. This collapses the former two-factor
 * posterior (MATreeLikelihood over sampled genotypes + LikelihoodReadCountModel) into a single
 * factor and removes the genotype StateNode and all its Gibbs operators.</p>
 *
 * <p>The read-count likelihoods are tiny, so each tip's partial vector is rescaled per-site by
 * its maximum before upload (max-normalisation); the summed log of those per-(tip,site) maxima
 * is held in {@link #globalLogOffset} and added back to {@code logP}.</p>
 *
 * <p>The {@code data} input is a constant scaffold alignment (e.g. {@code MutableAlignment} with
 * identity patterns) used only for taxa / site count / datatype / pattern indexing; its genotype
 * values are ignored. The read-count model referenced via {@link #readCountModelInput} must NOT
 * also be added to the posterior CompoundDistribution (it would be double-counted); it stays in
 * the model graph as this likelihood's input, so BEAST still drives its store/restore.</p>
 */
@Description("Tree likelihood that integrates single-cell genotypes out via read-count tip partials")
public class ReadCountTreeLikelihood extends TreeLikelihoodWithErrorFast {

    final public Input<LikelihoodReadCountModel> readCountModelInput = new Input<>(
            "readCountModel",
            "read-count model supplying per-genotype likelihoods P(reads | genotype) for the tip partials",
            Input.Validate.REQUIRED);

    private LikelihoodReadCountModel readCountModel;

    // alignment-taxon-index -> ReadCount-row-index
    private int[] alignToRC;

    // the 7 read-count parameters; a change in any forces a tip-partial rebuild
    private RealParameter[] rcParams;

    // summed log of the per-(tip,site) max-normalisation factors; added back to logP
    private double globalLogOffset = 0.0;
    private double storedGlobalLogOffset = 0.0;

    @Override
    public void initAndValidate() {
        readCountModel = readCountModelInput.get();
        // refresh cached Gamma tables before the first getLeafPartials traversal that
        // super.initAndValidate() triggers (the model's own calculateLogP never runs now)
        readCountModel.initialize();
        alignToRC = readCountModel.getAlignToRCIndex();
        rcParams = new RealParameter[]{
                readCountModel.epsilonInput.get(),
                readCountModel.deltaInput.get(),
                readCountModel.tInput.get(),
                readCountModel.vInput.get(),
                readCountModel.sInput.get(),
                readCountModel.w1Input.get(),
                readCountModel.w2Input.get(),
        };

        // super.initAndValidate() builds the initial tip partials via setPartials -> getLeafPartials,
        // accumulating globalLogOffset from 0.
        globalLogOffset = 0.0;
        super.initAndValidate();

        // sanity checks (see design review must-fixes #4, #5)
        Alignment data = dataInput.get();
        int stateCount = data.getDataType().getStateCount();
        if (stateCount != 10 && stateCount != 16) {
            throw new IllegalArgumentException(
                    "ReadCountTreeLikelihood expects a 10- or 16-state genotype datatype, got " + stateCount);
        }
        if (data.getPatternCount() != data.getSiteCount()) {
            throw new IllegalArgumentException(
                    "ReadCountTreeLikelihood requires an identity-pattern scaffold alignment "
                    + "(patternCount == siteCount); pattern compression would break the per-site "
                    + "read-count mapping. Use a MutableAlignment (or stripInvariant=false dense alignment).");
        }
        if (m_siteModel.getProportionInvariant() != 0.0) {
            throw new IllegalArgumentException(
                    "ReadCountTreeLikelihood does not support proportionInvariant > 0 "
                    + "(it would add an unscaled term to constant-pattern root slots). Set it to 0.");
        }
    }

    /**
     * Tip partial vector for one leaf: for each site (== pattern, identity), the per-genotype
     * read-count likelihood, max-normalised per site. Accumulates the per-site log-max into
     * {@link #globalLogOffset}.
     */
    @Override
    protected double[] getLeafPartials(Node node) {
        Alignment data = dataInput.get();
        int nrOfStates = data.getDataType().getStateCount();
        int nrOfPatterns = data.getPatternCount();
        double[] partials = new double[nrOfPatterns * nrOfStates];

        int alignIdx = getTaxonIndex(node.getID(), data);
        int rcIdx = alignToRC[alignIdx];

        double[] logRC = new double[nrOfStates];
        int i = 0;
        for (int p = 0; p < nrOfPatterns; p++) {
            int[] reads = readCountModel.getReadCount().getReadCounts(rcIdx, p);
            int coverage = readCountModel.getCoverage(alignIdx, p);
            double max = Double.NEGATIVE_INFINITY;
            for (int g = 0; g < nrOfStates; g++) {
                logRC[g] = readCountModel.logLiklihoodRC(g, reads, coverage, rcIdx);
                if (logRC[g] > max) max = logRC[g];
            }
            globalLogOffset += max;
            for (int g = 0; g < nrOfStates; g++) {
                partials[i++] = Math.exp(logRC[g] - max);
            }
        }
        return partials;
    }

    /**
     * Rebuild and re-upload every leaf's partials from the current read-count parameters,
     * recomputing globalLogOffset from scratch.
     */
    @Override
    public void updateLeafPartials() {
        globalLogOffset = 0.0;
        readCountModel.initialize();
        List<Node> leaves = treeInput.get().getExternalNodes();
        for (Node node : leaves) {
            int nodeId = node.getNr();
            likelihoodCore.setNodePartialsForUpdate(nodeId);
            BeerLikelihoodCore beer = (BeerLikelihoodCore) likelihoodCore;
            beer.setCurrentNodePartials(nodeId, getLeafPartials(node));
        }
    }

    @Override
    public double calculateLogP() {
        // super (Fast) re-runs updateLeafPartials() when the leaf-partial dirty flag is set,
        // which resets and re-accumulates globalLogOffset; otherwise globalLogOffset is unchanged.
        double core = super.calculateLogP();
        logP = core + globalLogOffset;
        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        // any read-count parameter change invalidates the tip partials -> full rebuild
        for (RealParameter p : rcParams) {
            if (p != null && p.somethingIsDirty()) {
                updateLeafPartials = true;
                hasDirt = Tree.IS_FILTHY;
                return true;
            }
        }
        // otherwise defer to the inherited tree / siteModel / branchRate / data checks
        return super.requiresRecalculation();
    }

    @Override
    public void store() {
        storedGlobalLogOffset = globalLogOffset;
        super.store();
    }

    @Override
    public void restore() {
        globalLogOffset = storedGlobalLogOffset;
        super.restore();
    }
}
