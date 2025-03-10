package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import mutablealignment.MATreeLikelihood;
import phylonco.lphy.evolution.readcountmodel.MutableAlignmentModel;

public class MutableAlignmentModelToBEAST implements GeneratorToBEAST<MutableAlignmentModel, MATreeLikelihood> {
    @Override
    public MATreeLikelihood generatorToBEAST(MutableAlignmentModel generator, BEASTInterface value, BEASTContext context) {
        MATreeLikelihood maTreeLikelihood = new MATreeLikelihood();
        return maTreeLikelihood;
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
