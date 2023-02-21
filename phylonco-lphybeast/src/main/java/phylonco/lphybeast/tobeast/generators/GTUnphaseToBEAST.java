package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import lphy.evolution.alignment.Alignment;
import lphy.graphicalModel.Generator;
import lphy.graphicalModel.GraphicalModelNode;
import lphy.graphicalModel.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.beast.evolution.likelihood.TreeLikelihoodWithError;
import phylonco.lphy.evolution.alignment.UnphaseGenotypeAlignment;

import java.util.List;

/**
 * This has to create TreeLikelihood.
 * A ~ PhyloCTMC(); E ~ ErrorModel(A); D = unphase(E);
 * @author Walter Xie
 */
public class GTUnphaseToBEAST implements GeneratorToBEAST<UnphaseGenotypeAlignment, TreeLikelihoodWithError>  {

    @Override
    public TreeLikelihoodWithError generatorToBEAST(UnphaseGenotypeAlignment generator, BEASTInterface value, BEASTContext context) {

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

        // TreeLikelihoodWithError should be already created for err alignment
        TreeLikelihoodWithError treeLikelihoodWithError =
                (TreeLikelihoodWithError) context.getBEASTObject(errAligGenerator);

        if (treeLikelihoodWithError == null)
            throw new IllegalArgumentException("Cannot find err alignment tree likelihood !");

        treeLikelihoodWithError.setInputValue("data", unphasedErrAlignment);

        treeLikelihoodWithError.initAndValidate();
        treeLikelihoodWithError.setID(unphasedErrAlignment.getID() + ".treeLikelihood");

        BEASTInterface errAlignment = context.getBEASTObject(errAlignmentInput);
        context.removeBEASTObject(errAlignment);
        // remove previous treeLikelihoodWithError added by GT16ErrorModelToBEAST
        context.removeBEASTObject(treeLikelihoodWithError);

        // logging
        context.addExtraLoggable(treeLikelihoodWithError);
        return treeLikelihoodWithError;
    }


    @Override
    public Class<UnphaseGenotypeAlignment> getGeneratorClass() {
        return UnphaseGenotypeAlignment.class;
    }

    @Override
    public Class<TreeLikelihoodWithError> getBEASTClass() {
        return TreeLikelihoodWithError.class;
    }
}
