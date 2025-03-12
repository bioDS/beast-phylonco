package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.sitemodel.SiteModel;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.likelihood.PhyloCTMC;
import lphy.core.model.GraphicalModelNode;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.generators.PhyloCTMCToBEAST;
import mutablealignment.MATreeLikelihood;
import phylonco.lphy.evolution.readcountmodel.MutableAlignmentModel;

import java.util.Objects;

public class MutableAlignmentModelToBEAST implements GeneratorToBEAST<MutableAlignmentModel, MATreeLikelihood> {
    @Override
    public MATreeLikelihood generatorToBEAST(MutableAlignmentModel generator, BEASTInterface value, BEASTContext context) {
//        assert value instanceof mutablealignment.MutableAlignment;
//        mutablealignment.MutableAlignment mutableAlignment = (mutablealignment.MutableAlignment) value;
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

        MATreeLikelihood maTreeLikelihood = new MATreeLikelihood();
        PhyloCTMCToBEAST.constructTreeAndBranchRate(phyloCTMC, maTreeLikelihood, context, true);
        SiteModel siteModel = PhyloCTMCToBEAST.constructSiteModel(phyloCTMC, context);
        maTreeLikelihood.setInputValue("data", value);
        maTreeLikelihood.setInputValue("siteModel", siteModel);
        maTreeLikelihood.initAndValidate();
        maTreeLikelihood.setID(value.getID() + ".treeLikelihood");

        // logging
        context.addExtraLoggable(maTreeLikelihood);

        removeOriginalTreeLikelihood(origAlignmentInput, phyloCTMC, context);


        return maTreeLikelihood;
    }




    private void removeOriginalTreeLikelihood(Value<Alignment> origAlignmentInput, PhyloCTMC phyloCTMC, BEASTContext context) {
        BEASTInterface beastOrigAlignment = context.getBEASTObject(origAlignmentInput);
        context.removeBEASTObject(beastOrigAlignment);

        BEASTInterface treeLikelihood = context.getBEASTObject(phyloCTMC);
        context.removeBEASTObject(treeLikelihood);
    }






    @Override
    public Class<MutableAlignmentModel> getGeneratorClass() {
        return MutableAlignmentModel.class;
    }

    @Override
    public Class<MATreeLikelihood> getBEASTClass() {
        return MATreeLikelihood.class;
    }
}
