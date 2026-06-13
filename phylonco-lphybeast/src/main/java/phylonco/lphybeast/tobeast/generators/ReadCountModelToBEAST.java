package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.spec.evolution.operator.AdaptableVarianceMultivariateNormalOperator;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Tree;
import beast.base.spec.inference.operator.Transform;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
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

import beast.base.spec.evolution.operator.UpDownOperator;
import beast.base.spec.type.Tensor;

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
 *
 * <p>beast3 note: LPhyBeast 2.0 dropped {@code BEASTContext.getAsRealParameter}; the read-count
 * parameters are fetched via {@link BEASTContext#getBEASTObject} (the classic {@code RealParameter}
 * produced for an LPhy value, exactly as {@code GT10ToBEAST} does for the nucleotide rates), so they
 * can drive the {@code Function}/{@code StateNode}-based AVMN and up/down operators below.</p>
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

        // beast3: getBEASTObject returns the spec param produced for each LPhy value
        // (RealScalarParam for the scalar distributions, RealVectorParam for s ~ ...replicates),
        // which match the LikelihoodReadCountModel spec-param inputs directly (no cast needed).
        BEASTInterface epsilonParam = context.getBEASTObject(epsilonValue);
        BEASTInterface deltaParam = context.getBEASTObject(deltaValue);
        BEASTInterface tParam = context.getBEASTObject(tValue);
        BEASTInterface vParam = context.getBEASTObject(vValue);
        BEASTInterface sParam = context.getBEASTObject(sValue);
        BEASTInterface w1Param = context.getBEASTObject(w1Value);
        BEASTInterface w2Param = context.getBEASTObject(w2Value);

        if (!(value instanceof ReadCount readCountData)) {
            throw new IllegalArgumentException("Require read count data");
        }

        // ---- build the read-count model as a per-genotype-likelihood HELPER (not a factor) ----
        // No genotype alignment: the helper derives its scaffold (cells, sites) from the read counts,
        // and the genotype datatype is set by ReadCountTreeLikelihood from the GT10/GT16 subst model.
        LikelihoodReadCountModel likelihoodReadCountModel = new LikelihoodReadCountModel();
        likelihoodReadCountModel.setInputValue("epsilon", epsilonParam);
        likelihoodReadCountModel.setInputValue("delta", deltaParam);
        likelihoodReadCountModel.setInputValue("t", tParam);
        likelihoodReadCountModel.setInputValue("v", vParam);
        likelihoodReadCountModel.setInputValue("s", sParam);
        likelihoodReadCountModel.setInputValue("w1", w1Param);
        likelihoodReadCountModel.setInputValue("w2", w2Param);
        likelihoodReadCountModel.setInputValue("readCount", readCountData);
        likelihoodReadCountModel.initAndValidate();

        // ---- build the integrated tree likelihood (the single posterior factor) ----
        // It builds its identity-pattern scaffold internally from the read counts (no data input).
        PhyloCTMC phyloCTMC = (PhyloCTMC) alignmentValue.getGenerator();

        ReadCountTreeLikelihood treeLikelihood = new ReadCountTreeLikelihood();
        PhyloCTMCToBEAST.constructTreeAndBranchRate(phyloCTMC, treeLikelihood, context, true);
        SiteModel siteModel = PhyloCTMCToBEAST.constructSiteModel(phyloCTMC, context);
        treeLikelihood.setInputValue("siteModel", siteModel);
        treeLikelihood.setInputValue("readCountModel", likelihoodReadCountModel);
        treeLikelihood.initAndValidate();
        treeLikelihood.setID((value.getID() != null ? value.getID() : "readCount") + ".treeLikelihood");
        // the scaffold was built internally during initAndValidate only to satisfy the parent
        // TreeLikelihood; drop it so the XML contains no genotype alignment (rebuilt from read counts)
        treeLikelihood.dropScaffoldInput();

        // ---- the genotype alignment A is integrated out: remove A and its tree likelihood from the
        //      model, so the XML contains no genotype alignment (the only data is the read counts) ----
        BEASTInterface maTreeLikelihood = context.getBEASTObject(phyloCTMC);
        if (maTreeLikelihood != null && maTreeLikelihood != treeLikelihood) {
            context.removeBEASTObject(maTreeLikelihood);
        }
        BEASTInterface alignmentParam = context.getBEASTObject(alignmentValue);
        if (alignmentParam == null) {
            alignmentParam = context.getBEASTObject(alignmentValue.getId());
        }
        if (alignmentParam != null) {
            context.removeBEASTObject(alignmentParam);
        }

        // ---- tree-aware genotype reconstruction logger (replaces the sampled-alignment log) ----
        addSampledGenotypeLogger(context, treeLikelihood, siteModel, likelihoodReadCountModel, readCountData);

        // beast3 TODO: the PR added custom AVMN (log-transform over t,v) and up/down operators for
        // the coverage parameters. Those use the classic Function/StateNode operator path, which the
        // spec params don't fit (RealScalarParam is not a beast.base.core.Function). For now rely on
        // LPhyBeast's default operator scheduling for the read-count parameters; re-add spec-package
        // operators here for better mixing once available.

        // add UpDownOperator to s (up), t (down) and v (down)
        List<Tensor> upParam = new ArrayList();
        upParam.add((RealVectorParam) sParam);
        List<Tensor> downParam = new ArrayList();
        downParam.add((RealScalarParam) tParam);
        downParam.add((RealScalarParam) vParam);
        addUpDownOperator(context, upParam, downParam, "readCountModel.stv.UpDownOperator");

        // add AVMN operator
//        List<Tensor> avmnParam = new ArrayList();
//        avmnParam.add((RealScalarParam) tParam);
//        avmnParam.add((RealVectorParam) sParam);
//        addAVMNOperator(context, avmnParam);

        return treeLikelihood;
    }

    // add up down operator
    private void addUpDownOperator(BEASTContext context, List<? extends Tensor> upArgs, List<? extends Tensor> downArgs, String id) {
        UpDownOperator operator = new UpDownOperator();
        operator.setID(id);
        for (Tensor arg: upArgs) {
            operator.setInputValue("up", upArgs);
        }
        for (Tensor arg: downArgs) {
            operator.setInputValue("down", upArgs);
        }
        operator.setInputValue("scaleFactor", 0.75);
        operator.setInputValue("weight", 10.0);
        context.addExtraOperator(operator);
    }

    // add AVMN operator
    // <operator id="AVMNOperator" spec="kernel.AdaptableVarianceMultivariateNormalOperator"
    // beta="0.05" burnin="100" initial="200" weight="15.0">
    // <transformations id="Transform$LogTransform" spec="operator.kernel.Transform$LogTransform">
    //  <f idref="t"/>
    //  <f idref="v"/>
    // </transformations>
    // </operator>
    //
    private void addAVMNOperator(BEASTContext context, List<Tensor> params) {
        AdaptableVarianceMultivariateNormalOperator operator = new AdaptableVarianceMultivariateNormalOperator();
        List<Transform> transforms = new ArrayList<>();
        for (Tensor p: params) {
            transforms.add(new Transform.LogTransform(p));
        }
        operator.setInputValue("transformations", transforms);
        operator.setInputValue("beta", 0.05);
        operator.setInputValue("burnin", 100);
        operator.setInputValue("initial", 200);
        operator.setInputValue("weight", 15.0);
        operator.setID("readCountModel.AVMN");
        context.addExtraOperator(operator);
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

    @Override
    public Class<ReadCountModel> getGeneratorClass() {
        return ReadCountModel.class;
    }

    @Override
    public Class<ReadCountTreeLikelihood> getBEASTClass() {
        return ReadCountTreeLikelihood.class;
    }
}
