package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;
import phylonco.lphy.evolution.readcountmodel.ReadCountModel;

public class ReadCountModelToBEAST implements GeneratorToBEAST<ReadCountModel, LikelihoodReadCountModel> {
    @Override
    public LikelihoodReadCountModel generatorToBEAST(ReadCountModel generator, BEASTInterface value, BEASTContext context) {
        LikelihoodReadCountModel likelihoodReadCountModel = new LikelihoodReadCountModel();
        //Get value from LPhy
        Value epsilonValue = generator.getParams().get("epsilon");
        Value deltaValue = (Value) generator.getParams().get("alpha").getGenerator().getParams().get("delta");
        Value tValue = (Value) generator.getParams().get("coverage").getGenerator().getParams().get("t");
        Value vValue = (Value) generator.getParams().get("coverage").getGenerator().getParams().get("v");
        Value sValue = (Value) generator.getParams().get("coverage").getGenerator().getParams().get("s");
        Value wValue = generator.getParams().get("w");
        Value alignmentValue = generator.getParams().get("D");

        BEASTInterface epsilonParam = context.getBEASTObject(epsilonValue);
        BEASTInterface deltaParam = context.getBEASTObject(deltaValue);
        BEASTInterface tParam = context.getBEASTObject(tValue);
        BEASTInterface vParam = context.getBEASTObject(vValue);
        BEASTInterface sParam = context.getBEASTObject(sValue);
        BEASTInterface wParam = context.getBEASTObject(wValue);
        BEASTInterface alignmentParam = context.getBEASTObject(alignmentValue);

        likelihoodReadCountModel.setInputValue("epsilon", epsilonParam);
        likelihoodReadCountModel.setInputValue("delta", deltaParam);
        likelihoodReadCountModel.setInputValue("t", tParam);
        likelihoodReadCountModel.setInputValue("v", vParam);
        likelihoodReadCountModel.setInputValue("s", sParam);
        likelihoodReadCountModel.setInputValue("w", wParam);
        likelihoodReadCountModel.setInputValue("alignment", alignmentParam);

        // beast readcount readCountData
        if (value instanceof ReadCount readCountData) {
            likelihoodReadCountModel.setInputValue("readCount", readCountData);
        } else {
            throw new IllegalArgumentException("Require read count data");
        }

        likelihoodReadCountModel.initAndValidate();

        return likelihoodReadCountModel;
    }

    @Override
    public Class<ReadCountModel> getGeneratorClass() {
        return ReadCountModel.class;
    }

    @Override
    public Class<LikelihoodReadCountModel> getBEASTClass() {
        return LikelihoodReadCountModel.class;
    }


}
