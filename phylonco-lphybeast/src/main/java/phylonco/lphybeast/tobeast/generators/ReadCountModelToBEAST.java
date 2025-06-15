package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.inference.StateNode;
import beast.base.inference.operator.UpDownOperator;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.likelihood.PhyloCTMC;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import mutablealignment.MATreeLikelihood;
import mutablealignment.MutableAlignment;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.beast.evolution.readcountmodel.GibbsAlignmentOperator;
import phylonco.beast.evolution.readcountmodel.GibbsSequenceOperator;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;
import mutablealignment.MutableAlignmentOperator;
import phylonco.lphy.evolution.readcountmodel.ReadCountModel;

import java.util.ArrayList;
import java.util.List;

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
        String w1ParamName = "w1";
        String w2ParamName = "w2";
        String alignmentParamName = "D";

        //get values from LPhy
        Value epsilonValue = generator.getParams().get(epsilonParamName);
        Value alphaValue = generator.getParams().get(alphaParamName);
        Value deltaValue = (Value) alphaValue.getGenerator().getParams().get(deltaParamName);
        Value coverageValue = generator.getParams().get(covParamName);
        Value tValue = (Value) coverageValue.getGenerator().getParams().get(tParamName);
        Value vValue = (Value) coverageValue.getGenerator().getParams().get(vParamName);
        Value sValue = (Value) coverageValue.getGenerator().getParams().get(sParamName);
        Value w1Value = generator.getParams().get(w1ParamName);
        Value w2Value = generator.getParams().get(w2ParamName);
        Value alignmentValue = generator.getParams().get(alignmentParamName);

        //convert LPhy values to Beast objects
        RealParameter epsilonParam = context.getAsRealParameter(epsilonValue);
        RealParameter deltaParam = context.getAsRealParameter(deltaValue);
        RealParameter tParam = context.getAsRealParameter(tValue);
        RealParameter vParam = context.getAsRealParameter(vValue);
        RealParameter sParam = context.getAsRealParameter(sValue);
        RealParameter w1Param = context.getAsRealParameter(w1Value);
        RealParameter w2Param = context.getAsRealParameter(w2Value);

        BEASTInterface alignmentParam = context.getBEASTObject(alignmentValue);
        if (alignmentParam == null) {
            alignmentParam = context.getBEASTObject(alignmentValue.getId());
        }

        //set input of likelihood model
        likelihoodReadCountModel.setInputValue("epsilon", epsilonParam);
        likelihoodReadCountModel.setInputValue("delta", deltaParam);
        likelihoodReadCountModel.setInputValue("t", tParam);
        likelihoodReadCountModel.setInputValue("v", vParam);
        likelihoodReadCountModel.setInputValue("s", sParam);
        likelihoodReadCountModel.setInputValue("w1", w1Param);
        likelihoodReadCountModel.setInputValue("w2", w2Param);
        likelihoodReadCountModel.setInputValue("alignment", alignmentParam);
        // beast readcount readCountData
        if (value instanceof ReadCount readCountData) {
            likelihoodReadCountModel.setInputValue("readCount", readCountData);
        } else {
            throw new IllegalArgumentException("Require read count data");
        }

        likelihoodReadCountModel.initAndValidate();

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

        context.addSkipOperator(sParam);
        context.addSkipLoggable(sParam);

        PhyloCTMC lphyTreeLikelihood = (PhyloCTMC) alignmentValue.getGenerator();
        MATreeLikelihood maTreeLikelihood = (MATreeLikelihood) context.getBEASTObject(lphyTreeLikelihood);

        //addGibbsAlignmentOperator(context, alignmentParam, maTreeLikelihood, likelihoodReadCountModel);
        addGibbsSequenceOperator(context, alignmentParam, maTreeLikelihood, likelihoodReadCountModel) ;
        context.addSkipOperator((StateNode) alignmentParam);

        return likelihoodReadCountModel;
    }



    private Transform addLogTransform(List<Function> functions) {

        Transform.LogTransform logTransform = new Transform.LogTransform();
        logTransform.setInputValue("f", functions);
        logTransform.initAndValidate();
        return logTransform;
    }

    private Transform addNoTransform(List<Function> functions) {

        Transform.NoTransform noTransform = new Transform.NoTransform();
        noTransform.setInputValue("f", functions);
        noTransform.initAndValidate();
        return noTransform;
    }

    private Transform addLogitTransform(List<Function> functions) {

        Transform.LogitTransform logitTransform = new Transform.LogitTransform();
        logitTransform.setInputValue("f", functions);
        logitTransform.initAndValidate();
        return logitTransform;
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

        // add operator
        context.addExtraOperator(operator);
        // skip default operator schedule

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

    private void addMutableAlignmentOperator(BEASTContext context, BEASTInterface alignment) {
        MutableAlignmentOperator operator = new MutableAlignmentOperator();
        operator.setInputValue("mutableAlignment", alignment);
        operator.initAndValidate();
        operator.setID("mutableAlignmentOperator");
        context.addExtraOperator(operator);
    }

    private void addGibbsSequenceOperator(BEASTContext context, BEASTInterface alignment, MATreeLikelihood maTreeLikelihood, LikelihoodReadCountModel likelihoodReadCountModel) {
        GibbsSequenceOperator operator = new GibbsSequenceOperator();
        MutableAlignment mutableAlignment = (MutableAlignment) alignment;
        operator.setInputValue("mutableAlignment", alignment);
        operator.setInputValue("maTreeLikelihood", maTreeLikelihood);
        operator.setInputValue("likelihoodReadCountModel", likelihoodReadCountModel);
        operator.setInputValue("weight", Math.pow(mutableAlignment.getTaxonCount(), 0.5));
        operator.initAndValidate();
        operator.setID("gibbsSequenceOperator");
        context.addExtraOperator(operator);
    }

    private void addGibbsAlignmentOperator(BEASTContext context, BEASTInterface alignment, MATreeLikelihood maTreeLikelihood, LikelihoodReadCountModel likelihoodReadCountModel) {
        GibbsAlignmentOperator operator = new GibbsAlignmentOperator();
        MutableAlignment mutableAlignment = (MutableAlignment) alignment;
        operator.setInputValue("mutableAlignment", alignment);
        operator.setInputValue("maTreeLikelihood", maTreeLikelihood);
        operator.setInputValue("likelihoodReadCountModel", likelihoodReadCountModel);
        operator.setInputValue("weight", Math.pow(1, 0.5));
        operator.initAndValidate();
        operator.setID("gibbsAlignmentOperator");
        context.addExtraOperator(operator);
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
