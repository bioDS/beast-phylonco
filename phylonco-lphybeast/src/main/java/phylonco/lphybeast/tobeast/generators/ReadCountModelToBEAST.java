package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.inference.parameter.RealParameter;
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
        String epsilonParamName = "epsilon";
        String alphaParamName = "alpha";
        String deltaParamName = "delta";
        String covParamName = "coverage";
        String tParamName = "t";
        String vParamName = "v";
        String sParamName = "s";
        String wParamName = "w";
        String alignmentParamName = "D";

        Value epsilonValue = generator.getParams().get(epsilonParamName);
        Value alphaValue = generator.getParams().get(alphaParamName);
        Value deltaValue = (Value) alphaValue.getGenerator().getParams().get(deltaParamName);
        Value coverageValue = generator.getParams().get(covParamName);
        Value tValue = (Value) coverageValue.getGenerator().getParams().get(tParamName);
        Value vValue = (Value) coverageValue.getGenerator().getParams().get(vParamName);
        Value sValue = (Value) coverageValue.getGenerator().getParams().get(sParamName);
        Value wValue = generator.getParams().get(wParamName);
        Value alignmentValue = generator.getParams().get(alignmentParamName);


        RealParameter epsilonParam = context.getAsRealParameter(epsilonValue);
        RealParameter deltaParam = context.getAsRealParameter(deltaValue);
        RealParameter tParam = context.getAsRealParameter(tValue);
        RealParameter vParam = context.getAsRealParameter(vValue);
        RealParameter sParam = context.getAsRealParameter(sValue);
        RealParameter wParam = context.getAsRealParameter(wValue);
        BEASTInterface alignmentParam = context.getBEASTObject(alignmentValue);
        if (alignmentParam == null) {
            alignmentParam = context.getBEASTObject(alignmentValue.getId());
        }

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
