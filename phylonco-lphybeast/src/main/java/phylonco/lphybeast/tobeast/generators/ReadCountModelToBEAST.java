package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.inference.StateNode;
import beast.base.inference.operator.UpDownOperator;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.likelihood.PhyloCTMC;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.generators.PhyloCTMCToBEAST;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.beast.evolution.likelihood.ReadCountTreeLikelihood;
import phylonco.beast.evolution.likelihood.SampledGenotypeLogger;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;
import phylonco.lphy.evolution.readcountmodel.ReadCountModel;

import java.util.ArrayList;
import java.util.List;

/**
 * Translates an LPhy {@code r ~ ReadCountModel(D=A, ...)} into the integrated-genotype BEAST model:
 * the latent genotype alignment {@code A} is integrated out analytically via read-count tip partials
 * rather than sampled. The single posterior factor is a {@link ReadCountTreeLikelihood} (the genotype
 * tree likelihood with read-count tip partials); the {@link LikelihoodReadCountModel} is wired in only
 * as its per-genotype-likelihood helper (NOT a separate factor). The original
 * {@code mutablealignment.MATreeLikelihood} is removed and no genotype Gibbs operators are created —
 * {@code A} stays as a constant identity-pattern scaffold (skipped from operators and logging).
 *
 * <p>Mirrors the teardown-and-replace pattern of {@link GT16ErrorModelToBEAST}.</p>
 */
public class ReadCountModelToBEAST implements GeneratorToBEAST<ReadCountModel, ReadCountTreeLikelihood> {
    @Override
    public ReadCountTreeLikelihood generatorToBEAST(ReadCountModel generator, BEASTInterface value, BEASTContext context) {
        // ---- pull the LPhy parameter values ----
        Value epsilonValue = generator.getParams().get("epsilon");
        Value alphaValue = generator.getParams().get("alpha");
        Value deltaValue = (Value) alphaValue.getGenerator().getParams().get("delta");
        Value coverageValue = generator.getParams().get("coverage");
        Value tValue = (Value) coverageValue.getGenerator().getParams().get("t");
        Value vValue = (Value) coverageValue.getGenerator().getParams().get("v");
        Value sValue = (Value) coverageValue.getGenerator().getParams().get("s");
        Value w1Value = generator.getParams().get("w1");
        Value w2Value = generator.getParams().get("w2");
        Value alignmentValue = generator.getParams().get("D");

        RealParameter epsilonParam = context.getAsRealParameter(epsilonValue);
        RealParameter deltaParam = context.getAsRealParameter(deltaValue);
        RealParameter tParam = context.getAsRealParameter(tValue);
        RealParameter vParam = context.getAsRealParameter(vValue);
        RealParameter sParam = context.getAsRealParameter(sValue);
        RealParameter w1Param = context.getAsRealParameter(w1Value);
        RealParameter w2Param = context.getAsRealParameter(w2Value);

        // the genotype alignment A: kept as a constant identity-pattern scaffold (values ignored)
        BEASTInterface alignmentParam = context.getBEASTObject(alignmentValue);
        if (alignmentParam == null) {
            alignmentParam = context.getBEASTObject(alignmentValue.getId());
        }

        if (!(value instanceof ReadCount readCountData)) {
            throw new IllegalArgumentException("Require read count data");
        }

        // ---- build the read-count model as a per-genotype-likelihood HELPER (not a factor) ----
        LikelihoodReadCountModel likelihoodReadCountModel = new LikelihoodReadCountModel();
        likelihoodReadCountModel.setInputValue("epsilon", epsilonParam);
        likelihoodReadCountModel.setInputValue("delta", deltaParam);
        likelihoodReadCountModel.setInputValue("t", tParam);
        likelihoodReadCountModel.setInputValue("v", vParam);
        likelihoodReadCountModel.setInputValue("s", sParam);
        likelihoodReadCountModel.setInputValue("w1", w1Param);
        likelihoodReadCountModel.setInputValue("w2", w2Param);
        likelihoodReadCountModel.setInputValue("alignment", alignmentParam);
        likelihoodReadCountModel.setInputValue("readCount", readCountData);
        likelihoodReadCountModel.initAndValidate();

        // ---- build the integrated tree likelihood (the single posterior factor) ----
        PhyloCTMC phyloCTMC = (PhyloCTMC) alignmentValue.getGenerator();

        ReadCountTreeLikelihood treeLikelihood = new ReadCountTreeLikelihood();
        treeLikelihood.setInputValue("data", alignmentParam);
        PhyloCTMCToBEAST.constructTreeAndBranchRate(phyloCTMC, treeLikelihood, context, true);
        SiteModel siteModel = PhyloCTMCToBEAST.constructSiteModel(phyloCTMC, context);
        treeLikelihood.setInputValue("siteModel", siteModel);
        treeLikelihood.setInputValue("readCountModel", likelihoodReadCountModel);
        treeLikelihood.initAndValidate();
        treeLikelihood.setID(alignmentParam.getID() + ".treeLikelihood");

        // ---- remove the original MATreeLikelihood; A is no longer sampled ----
        BEASTInterface maTreeLikelihood = context.getBEASTObject(phyloCTMC);
        if (maTreeLikelihood != null && maTreeLikelihood != treeLikelihood) {
            context.removeBEASTObject(maTreeLikelihood);
        }
        // keep A as a constant scaffold: no operators, not logged
        context.addSkipOperator((StateNode) alignmentParam);
        context.addSkipLoggable((StateNode) alignmentParam);

        // ---- tree-aware genotype reconstruction logger (replaces the sampled-alignment log) ----
        addSampledGenotypeLogger(context, treeLikelihood, siteModel, likelihoodReadCountModel, readCountData);

        // ---- operators for the read-count / coverage parameters ----
        List<Transform> transforms = new ArrayList<>();
        List<Function> logFunctions = new ArrayList<>();
        logFunctions.add(tParam);
        logFunctions.add(vParam);
        transforms.add(addLogTransform(logFunctions));
        addAVMNOperator(context, transforms);

        List<StateNode> upStates = new ArrayList<>();
        List<StateNode> downStates = new ArrayList<>();
        downStates.add(tParam);
        downStates.add(vParam);
        addUpDownOperator(context, upStates, downStates);

        return treeLikelihood;
    }

    private void addSampledGenotypeLogger(BEASTContext context, ReadCountTreeLikelihood treeLikelihood,
                                          SiteModel siteModel, LikelihoodReadCountModel readCountModel,
                                          ReadCount readCount) {
        SampledGenotypeLogger logger = new SampledGenotypeLogger();
        logger.setInputValue("tree", (Tree) treeLikelihood.treeInput.get());
        logger.setInputValue("siteModel", siteModel);
        BranchRateModel branchRateModel = treeLikelihood.branchRateModelInput.get();
        if (branchRateModel != null) {
            logger.setInputValue("branchRateModel", branchRateModel);
        }
        logger.setInputValue("readCountModel", readCountModel);
        logger.setInputValue("readCount", readCount);
        logger.initAndValidate();
        logger.setID("sampledGenotype");
        context.addExtraLoggable(logger);
    }

    private Transform addLogTransform(List<Function> functions) {
        Transform.LogTransform logTransform = new Transform.LogTransform();
        logTransform.setInputValue("f", functions);
        logTransform.initAndValidate();
        return logTransform;
    }

    private void addAVMNOperator(BEASTContext context, List<Transform> transforms) {
        AdaptableVarianceMultivariateNormalOperator operator = new AdaptableVarianceMultivariateNormalOperator();
        operator.setInputValue("weight", 8.0);
        operator.setInputValue("coefficient", 1.0);
        operator.setInputValue("scaleFactor", 1.0);
        operator.setInputValue("beta", 0.05);
        operator.setInputValue("initial", 200 * transforms.size());
        operator.setInputValue("burnin", 100 * transforms.size());
        operator.setInputValue("every", 1);
        operator.setInputValue("allowNonsense", false);
        operator.setInputValue("transformations", transforms);
        operator.initAndValidate();
        operator.setID("AVMN");
        context.addExtraOperator(operator);
    }

    private void addUpDownOperator(BEASTContext context, List<StateNode> upStateNode, List<StateNode> downStateNode) {
        UpDownOperator operator = new UpDownOperator();
        operator.setInputValue("up", upStateNode);
        operator.setInputValue("down", downStateNode);
        operator.setInputValue("weight", 3.0);
        operator.setInputValue("scaleFactor", 0.75);
        operator.initAndValidate();
        operator.setID("upDownOperator");
        context.addExtraOperator(operator);
    }

    @Override
    public Class<ReadCountModel> getGeneratorClass() {
        return ReadCountModel.class;
    }

    @Override
    public Class<ReadCountTreeLikelihood> getBEASTClass() {
        return ReadCountTreeLikelihood.class;
    }
}
