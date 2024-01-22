package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.likelihood.PhyloCTMC;
import lphy.core.model.GraphicalModelNode;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.generators.PhyloCTMCToBEAST;
import phylonco.beast.evolution.errormodel.ErrorModel;
import phylonco.beast.evolution.likelihood.TreeLikelihoodWithError;
import phylonco.beast.evolution.likelihood.TreeLikelihoodWithErrorFast;
import phylonco.lphy.evolution.alignment.GT16ErrorModel;

import java.util.Objects;

/**
 * This has to create TreeLikelihood
 * @author Walter Xie
 */
public class GT16ErrorModelToBEAST implements GeneratorToBEAST<GT16ErrorModel, TreeLikelihoodWithError> {

    @Override
    public TreeLikelihoodWithError generatorToBEAST(GT16ErrorModel generator, BEASTInterface value, BEASTContext context) {

        assert value instanceof beast.base.evolution.alignment.Alignment;
        beast.base.evolution.alignment.Alignment errAlignment = (beast.base.evolution.alignment.Alignment) value;

        phylonco.beast.evolution.errormodel.GT16ErrorModel gt16ErrorModel =
                new phylonco.beast.evolution.errormodel.GT16ErrorModel();

        DataType beastDataType = errAlignment.getDataType();
        gt16ErrorModel.setInputValue("datatype", beastDataType);

        RealParameter deltaParam = context.getAsRealParameter(generator.getDelta());
        gt16ErrorModel.setInputValue("delta", deltaParam);
        RealParameter epsilonParam = context.getAsRealParameter(generator.getEpsilon());
        gt16ErrorModel.setInputValue("epsilon", epsilonParam);
        gt16ErrorModel.initAndValidate();

        // Temporary solution to rm parent alignment if there is a child alignment created from it,
        // e.g. original alignment creates err alignment
        // A ~ PhyloCTMC(); E ~ ErrorModel(A);
        PhyloCTMC phyloCTMC = null;
        Value<Alignment> origAlignmentInput = null;
        for (GraphicalModelNode<?> input : Objects.requireNonNull(generator.getInputs())) {
            if (input instanceof Value && input.value() instanceof Alignment) {
                origAlignmentInput = (Value<Alignment>) input;
                phyloCTMC = (PhyloCTMC) origAlignmentInput.getGenerator();
                break;
            }
        }

        if (phyloCTMC == null) {
            throw new IllegalArgumentException("Cannot find err alignment and PhyloCTMC !");
        }

        TreeLikelihoodWithError treeLikelihoodWithError =
                getTreeLikelihoodWithError(errAlignment, gt16ErrorModel, phyloCTMC, context);

        // logging
        context.addExtraLoggable(treeLikelihoodWithError);

        removeOriginalTreeLikelihood(origAlignmentInput, phyloCTMC, context);

        return treeLikelihoodWithError;
    }

    private void removeOriginalTreeLikelihood(Value<Alignment> origAlignmentInput, PhyloCTMC phyloCTMC, BEASTContext context) {
        BEASTInterface beastOrigAlignment = context.getBEASTObject(origAlignmentInput);
        context.removeBEASTObject(beastOrigAlignment);

        BEASTInterface treeLikelihood = context.getBEASTObject(phyloCTMC);
        context.removeBEASTObject(treeLikelihood);
    }


    private TreeLikelihoodWithError getTreeLikelihoodWithError(beast.base.evolution.alignment.Alignment errAlignment,
                                                               ErrorModel errorModel, PhyloCTMC phyloCTMC, BEASTContext context) {
        TreeLikelihoodWithErrorFast treeLikelihoodWithError = new TreeLikelihoodWithErrorFast();

        treeLikelihoodWithError.setInputValue("data", errAlignment);

        // branch rate operators already created by generic TreeLikeihood
        PhyloCTMCToBEAST.constructTreeAndBranchRate(phyloCTMC, treeLikelihoodWithError, context, true);

        SiteModel siteModel = PhyloCTMCToBEAST.constructSiteModel(phyloCTMC, context);

        treeLikelihoodWithError.setInputValue("siteModel", siteModel);
        treeLikelihoodWithError.setInputValue("errorModel", errorModel);
        treeLikelihoodWithError.setInputValue("useTipsEmpirical", false);

        treeLikelihoodWithError.initAndValidate();
        treeLikelihoodWithError.setID(errAlignment.getID() + ".treeLikelihood");

        return treeLikelihoodWithError;
    }


    @Override
    public Class<GT16ErrorModel> getGeneratorClass() {
        return GT16ErrorModel.class;
    }

    @Override
    public Class<TreeLikelihoodWithError> getBEASTClass() {
        return TreeLikelihoodWithError.class;
    }
}
