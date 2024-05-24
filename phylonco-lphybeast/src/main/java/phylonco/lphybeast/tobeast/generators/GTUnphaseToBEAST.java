package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.ThreadedTreeLikelihood;
import lphy.base.evolution.alignment.Alignment;
import lphy.core.model.Generator;
import lphy.core.model.GraphicalModelNode;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.beast.evolution.likelihood.TreeLikelihoodWithError;
import phylonco.lphy.evolution.alignment.UnphaseGenotypeAlignment;

import java.util.List;

/**
 * This has to create TreeLikelihood.
 * A ~ PhyloCTMC();
 * E ~ ErrorModel(A);
 * D = unphase(E);
 * @author Walter Xie
 * @author Kylie Chen
 * @author Yuan Xu
 */
public class GTUnphaseToBEAST implements GeneratorToBEAST<UnphaseGenotypeAlignment, GenericTreeLikelihood>  {

    @Override
    public GenericTreeLikelihood generatorToBEAST(UnphaseGenotypeAlignment generator, BEASTInterface value, BEASTContext context) {

        assert value instanceof beast.base.evolution.alignment.Alignment;
        beast.base.evolution.alignment.Alignment unphasedErrAlignment = (beast.base.evolution.alignment.Alignment)value;

        Value<Alignment> errAlignmentInput = null;
        Generator<?>  errAligGenerator = null;
        List<GraphicalModelNode<?>> inputs = generator.getInputs();
        for (GraphicalModelNode<?> input : inputs) {
            if (input instanceof Value && input.value() instanceof Alignment) {
                errAlignmentInput = (Value<Alignment>) input;
                errAligGenerator = errAlignmentInput.getGenerator();
                break;
            }
        }

        GenericTreeLikelihood treeLikelihood = null;

        // only cast if TreeLikelihoodWithError if using an error model
        if (context.getBEASTObject(errAligGenerator) instanceof TreeLikelihoodWithError) {
            treeLikelihood = (TreeLikelihoodWithError) context.getBEASTObject(errAligGenerator);
        } else {
            treeLikelihood = (ThreadedTreeLikelihood) context.getBEASTObject(errAligGenerator);
        }

        if (treeLikelihood == null)
            throw new IllegalArgumentException("Cannot find err alignment tree likelihood !");

        treeLikelihood.setInputValue("data", unphasedErrAlignment);



        //Add this open the ambiguities function
        treeLikelihood.setInputValue("useAmbiguities", true);
//        treeLikelihood.setInputValue("useTipLikelihoods", true);


        treeLikelihood.initAndValidate();
        treeLikelihood.setID(unphasedErrAlignment.getID() + ".treeLikelihood");

        BEASTInterface errAlignment = context.getBEASTObject(errAlignmentInput);
        context.removeBEASTObject(errAlignment);
        // remove previous treeLikelihood added by GT16ErrorModelToBEAST
        context.removeBEASTObject(treeLikelihood);

        // logging
        context.addExtraLoggable(treeLikelihood);
        return treeLikelihood;
    }


    @Override
    public Class<UnphaseGenotypeAlignment> getGeneratorClass() {
        return UnphaseGenotypeAlignment.class;
    }

    @Override
    public Class<GenericTreeLikelihood> getBEASTClass() {
        return GenericTreeLikelihood.class;
    }
}
