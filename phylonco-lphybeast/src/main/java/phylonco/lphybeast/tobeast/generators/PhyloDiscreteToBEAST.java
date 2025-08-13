package phylonco.lphybeast.tobeast.generators;

import NestedBD.evolution.likelihood.DiploidOriginLikelihood;
import beast.base.core.BEASTInterface;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.RealParameter;
import lphy.base.distribution.LogNormal;
import lphy.base.distribution.UCLNMean1;
import lphy.base.evolution.branchrate.LocalBranchRates;
import lphy.base.evolution.branchrate.LocalClock;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.Generator;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.vectorization.IID;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.operators.DefaultOperatorStrategy;
import phylonco.lphy.evolution.copynumbermodel.MarkovTraitEvolution;
import phylonco.lphy.evolution.copynumbermodel.PhyloDiscrete;
import phylonco.lphybeast.tobeast.values.CopyNumberBDToBEAST;

import static lphybeast.tobeast.generators.PhyloCTMCToBEAST.addORCOperators;
import static lphybeast.tobeast.generators.PhyloCTMCToBEAST.getClockRateParam;

public class PhyloDiscreteToBEAST implements GeneratorToBEAST<PhyloDiscrete, DiploidOriginLikelihood> {

    @Override
    public DiploidOriginLikelihood generatorToBEAST(PhyloDiscrete generator, BEASTInterface value, BEASTContext context) {
        // Create the DiploidOriginLikelihood
        DiploidOriginLikelihood likelihood = new DiploidOriginLikelihood();

        // Get the beast integer alignment
        if (value instanceof beast.base.evolution.alignment.Alignment data) {
            // TODO: check data is type integer alignment
            likelihood.setInputValue("data", data); // Set alignment data
        } else {
            throw new IllegalArgumentException("Require alignment data");
        }
        //
        constructTreeAndBranchRate(generator, likelihood, context);

        // Set origin time
        likelihood.setInputValue("origtime", new RealParameter("0.0"));

        // Use the default nstates
        int nstates = CopyNumberBDToBEAST.DEFAULT_NSTATES;
        likelihood.setInputValue("nstates", new RealParameter(String.valueOf(nstates)));

        // Get the model
        Value<MarkovTraitEvolution<Integer>> modelValue = generator.getModel();
        SubstitutionModel evolutionModel = (SubstitutionModel) context.getBEASTObject(modelValue);
        // Create a site model with the substitution model
        beast.base.evolution.sitemodel.SiteModel siteModel = new beast.base.evolution.sitemodel.SiteModel();
        siteModel.setInputValue("substModel", evolutionModel); // Set substModel
        siteModel.initAndValidate();
        // Set siteModel
        likelihood.setInputValue("siteModel", siteModel);
        likelihood.initAndValidate();

        return likelihood;
    }

        // Create tree and clock rate inside this tree likelihood.
    public static void constructTreeAndBranchRate(PhyloDiscrete phyloDiscrete, DiploidOriginLikelihood likelihood, BEASTContext context) {
        constructTreeAndBranchRate(phyloDiscrete, likelihood, context, false);
    }

    public static void constructTreeAndBranchRate (PhyloDiscrete PhyloDiscrete, DiploidOriginLikelihood
            likelihood, BEASTContext context,boolean skipBranchOperators){
        // Get and Set the tree
        Value<TimeTree> treeValue = PhyloDiscrete.getTree();
        Tree tree = (Tree) context.getBEASTObject(treeValue);
        likelihood.setInputValue("tree", tree);

        // Handle clock rate
        Value<Number> clockRateValue = PhyloDiscrete.getClockRate();
        // Get clock rate
        RealParameter clockRateParam = getClockRateParam(clockRateValue, context);
        // updown op when estimating clock.rate
        if (clockRateValue instanceof RandomVariable && treeValue instanceof RandomVariable && skipBranchOperators == false) {
            DefaultOperatorStrategy.addUpDownOperator(tree, clockRateParam, context);
        }
        // relaxed or local clock
        Value<Double[]> branchRates = PhyloDiscrete.getBranchRates();
        if (branchRates != null) {
            /**
             * 1. Use ORC package (relaxed clock) where UCLNMean1 is used in LPhy
             * 2. Keep the alternative option to use IID LogNormal on branch rates in LPhy, but XML is not recommended.
             */
            Generator generator = branchRates.getGenerator();
            if (generator instanceof UCLNMean1 ucln) {

                // UCLNRelaxedClockToBEAST: the mean of log-normal distr on branch rates in real space is fixed to 1.
                UCRelaxedClockModel relaxedClockModel = (UCRelaxedClockModel) context.getBEASTObject(generator);
                likelihood.setInputValue("branchRateModel", relaxedClockModel);

                if (skipBranchOperators == false) {
                    addORCOperators(tree, relaxedClockModel, context);
                }

            } else if (generator instanceof IID &&
                    ((IID<?>) generator).getBaseDistribution() instanceof LogNormal) {
                LoggerUtils.log.warning("To use ORC package, please use UCLN_Mean1 in your lphy script !");

                // simpleRelaxedClock.lphy
                UCRelaxedClockModel relaxedClockModel = new UCRelaxedClockModel();

                Prior logNormalPrior = (Prior) context.getBEASTObject(generator);

                RealParameter beastBranchRates = context.getAsRealParameter(branchRates);

                relaxedClockModel.setInputValue("rates", beastBranchRates);
                relaxedClockModel.setInputValue("tree", tree);
                relaxedClockModel.setInputValue("distr", logNormalPrior.distInput.get());
                relaxedClockModel.setID(branchRates.getCanonicalId() + ".model");
                relaxedClockModel.initAndValidate();
                likelihood.setInputValue("branchRateModel", relaxedClockModel);

                if (skipBranchOperators == false) {
                    addORCOperators(tree, relaxedClockModel, context);
                }

            } else if (generator instanceof LocalBranchRates) {
                likelihood.setInputValue("branchRateModel", context.getBEASTObject(generator));
            } else if (generator instanceof LocalClock) {
                likelihood.setInputValue("branchRateModel", context.getBEASTObject(generator));
            } else {
                throw new UnsupportedOperationException("Only localBranchRates and lognormally distributed branchRates currently supported for LPhyBEAST !");
            }

        } else {
            StrictClockModel clockModel = new StrictClockModel();
            clockModel.setInputValue("clock.rate", clockRateParam);
            likelihood.setInputValue("branchRateModel", clockModel);
        }
    }

    @Override
    public Class<PhyloDiscrete> getGeneratorClass() {
        return PhyloDiscrete.class;
    }

    @Override
    public Class<DiploidOriginLikelihood> getBEASTClass() {
        return DiploidOriginLikelihood.class;
    }
}
