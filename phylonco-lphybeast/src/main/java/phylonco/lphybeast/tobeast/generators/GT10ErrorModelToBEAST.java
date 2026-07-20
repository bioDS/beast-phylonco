package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.datatype.DataType;

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
import phylonco.lphy.evolution.alignment.GT10ErrorModel;

import beast.base.spec.evolution.sitemodel.SiteModel;

import java.util.Objects;


public class GT10ErrorModelToBEAST implements GeneratorToBEAST<GT10ErrorModel, TreeLikelihoodWithError> {

    @Override
    public TreeLikelihoodWithError generatorToBEAST(GT10ErrorModel lphyErrorModel, BEASTInterface value, BEASTContext context) {

        assert value instanceof beast.base.evolution.alignment.Alignment;
        beast.base.evolution.alignment.Alignment errAlignment = (beast.base.evolution.alignment.Alignment) value;

        phylonco.beast.evolution.errormodel.GT10ErrorModel beastGT10ErrorModel =
                new phylonco.beast.evolution.errormodel.GT10ErrorModel();

        DataType beastDataType = errAlignment.getDataType();
        beastGT10ErrorModel.setInputValue("datatype", beastDataType);

        beastGT10ErrorModel.setInputValue("delta", context.getBEASTObject(lphyErrorModel.getDelta()));
        beastGT10ErrorModel.setInputValue("epsilon", context.getBEASTObject(lphyErrorModel.getEpsilon()));
        beastGT10ErrorModel.initAndValidate();

        // Temporary solution to rm parent alignment if there is a child alignment created from it,
        // e.g. original alignment creates err alignment
        // A ~ PhyloCTMC(); E ~ ErrorModel(A);
        PhyloCTMC phyloCTMC = null;
        Value<Alignment> origAlignmentInput = null;
        for (GraphicalModelNode<?> input : Objects.requireNonNull(lphyErrorModel.getInputs())) {
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
                getTreeLikelihoodWithError(errAlignment, beastGT10ErrorModel, phyloCTMC, context);

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
    public Class<GT10ErrorModel> getGeneratorClass() {
        return GT10ErrorModel.class;
    }

    @Override
    public Class<TreeLikelihoodWithError> getBEASTClass() {
        return TreeLikelihoodWithError.class;
    }
}
