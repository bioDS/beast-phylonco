package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.evolution.operator.AdaptableOperatorSampler;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;

import beast.base.inference.operator.DeltaExchangeOperator;
import beast.base.inference.operator.SwapOperator;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.lphy.evolution.substitutionmodel.GT10;

import beast.base.spec.type.Simplex;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.evolution.substitutionmodel.Frequencies;

import java.util.ArrayList;
import java.util.List;

/**
 * This has to create TreeLikelihood.
 * A ~ PhyloCTMC();
 * D ~ ErrorModel(A);
 */
public class GT10ToBEAST implements GeneratorToBEAST<GT10, phylonco.beast.evolution.substitutionmodel.GT10> {

    @Override
    public phylonco.beast.evolution.substitutionmodel.GT10 generatorToBEAST(GT10 gt10, BEASTInterface value, BEASTContext context) {

        phylonco.beast.evolution.substitutionmodel.GT10 beastGT10 = new phylonco.beast.evolution.substitutionmodel.GT10();

        Value<Double[]> rates = gt10.getRates();
        Value<Double[]> freqs = gt10.getFreq();

        // beast3: getBEASTObject returns the spec param produced for the LPhy value
        // (a Dirichlet rates vector -> SimplexParam, which is a RealVectorParam).
        RealVectorParam ratesParameter = (RealVectorParam) context.getBEASTObject(rates);
        ratesParameter.setInputValue("keys", "AC AG AT CG CT GT");
        ratesParameter.initAndValidate();

        Frequencies freqsParameter = new Frequencies((Simplex) context.getBEASTObject(gt10.getFreq()));

        beastGT10.setInputValue("nucRates", ratesParameter);
        beastGT10.setInputValue("frequencies", freqsParameter);

        beastGT10.initAndValidate();

        // beast3 TODO: the custom AVMN operators on the spec Simplex rates/frequencies need the
        // spec-package operators (the classic AVMN/Transform Function path can't take a Simplex).
        // Until then, rely on LPhyBeast's default operator scheduling for these parameters.

        return beastGT10;
    }

    private void addDeltaExchangeOperator(RealParameter parameter, BEASTContext context) {
        DeltaExchangeOperator operator = new DeltaExchangeOperator();
        operator.setInputValue("parameter", parameter);
        operator.setInputValue("weight", 2*context.getOperatorWeight(parameter.getDimension() - 1));
        operator.setInputValue("autoOptimize", true);
        operator.initAndValidate();
        operator.setID(parameter.getID() + ".deltaExchange");
        // add operator
        context.addExtraOperator(operator);
        // skip default operator schedule
        context.addSkipOperator(parameter);
    }

    private void addExtraDeltaExchangeOperators(Frequencies freqsParameter, BEASTContext context) {
        // heterozygous pairs
        String[] filterArray = {
                "0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0",
                "0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0",
                "0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0",
                "0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0",
                "0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0",
                "0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0"
        };

        IntegerParameter[] pairs = new IntegerParameter[filterArray.length];
        for (int i = 0; i < filterArray.length; i++) {
            String filter = filterArray[i];
            pairs[i] = new IntegerParameter(filter);
        };

        for (int i = 0; i < pairs.length; i++) {
            DeltaExchangeOperator operator = new DeltaExchangeOperator();
            operator.setInputValue("parameter", freqsParameter.frequenciesInput.get());
            operator.setInputValue("autoOptimize", true);
            operator.setInputValue("delta", 1.0 / 16);
            operator.setInputValue("weightvector", pairs[i]);
            operator.setInputValue("weight", "2.0"); // set operator weight to 1
            operator.setID("deltaExchangePair" + (i + 1));
            operator.initAndValidate();
            // add operator
            context.addExtraOperator(operator);
        }
    }

    private void addExtraSwapOperators(Frequencies freqsParameter, BEASTContext context) {
        // heterozygous pairs
        String[] filterArray = {
                "0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0",
                "0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0",
                "0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0",
                "0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0",
                "0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0",
                "0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0"
        };

        BooleanParameter[] pairs = new BooleanParameter[filterArray.length];
        for (int i = 0; i < filterArray.length; i++) {
            String filter = filterArray[i];
            pairs[i] = new BooleanParameter(filter);
        };

        for (int i = 0; i < pairs.length; i++) {
            SwapOperator operator = new SwapOperator();
            operator.setInputValue("parameter", freqsParameter.frequenciesInput.get());
            operator.setInputValue("filter", pairs[i]);
            operator.setInputValue("weight", "0.5"); // set operator weight to 0.5
            operator.setID("swapPair" + (i + 1));
            operator.initAndValidate();
            // add operator
            context.addExtraOperator(operator);
        }
    }

    private void addAVMNOperator(RealParameter parameter, BEASTContext context) {
        AdaptableVarianceMultivariateNormalOperator adaptableVarianceMultivariateNormalOperator = new AdaptableVarianceMultivariateNormalOperator();
        List<Transform> transforms = new ArrayList<>();
        List<Function> logFunctions = new ArrayList<>();
        logFunctions.add(parameter);
        Transform.LogConstrainedSumTransform logConstrainedSumTransform = new Transform.LogConstrainedSumTransform();
        logConstrainedSumTransform.setInputValue("f", logFunctions);
        logConstrainedSumTransform.initAndValidate();
        transforms.add(logConstrainedSumTransform);
        adaptableVarianceMultivariateNormalOperator.setInputValue("weight", 3*context.getOperatorWeight(parameter.getDimension() - 1));
        adaptableVarianceMultivariateNormalOperator.setInputValue("coefficient", 1.0);
        adaptableVarianceMultivariateNormalOperator.setInputValue("scaleFactor", 1.0);
        adaptableVarianceMultivariateNormalOperator.setInputValue("beta", 0.05);
        adaptableVarianceMultivariateNormalOperator.setInputValue("initial", 200 * transforms.size());
        adaptableVarianceMultivariateNormalOperator.setInputValue("burnin", 100 * transforms.size());
        adaptableVarianceMultivariateNormalOperator.setInputValue("every", 1);
        adaptableVarianceMultivariateNormalOperator.setInputValue("allowNonsense", false);
        adaptableVarianceMultivariateNormalOperator.setInputValue("transformations", transforms);
        adaptableVarianceMultivariateNormalOperator.initAndValidate();
        adaptableVarianceMultivariateNormalOperator.setID(parameter.getID() + ".AVMN");
        context.addExtraOperator(adaptableVarianceMultivariateNormalOperator);
        context.addSkipOperator(parameter);
    }


    private void addAdaptableOperatorSampler(RealParameter parameter, BEASTContext context) {
        double weight = 0.0;
        // skip default operator schedule
        context.addSkipOperator(parameter);
        //add AdaptableOperatorSampler
        AdaptableOperatorSampler adaptableOperatorSampler = new AdaptableOperatorSampler();
        adaptableOperatorSampler.setInputValue("parameter", parameter);

        //add DeltaExchangeOperator
        DeltaExchangeOperator deltaExchangeOperator = new DeltaExchangeOperator();
        deltaExchangeOperator.setInputValue("parameter", parameter);
        deltaExchangeOperator.setInputValue("weight", context.getOperatorWeight(parameter.getDimension() - 1));
        deltaExchangeOperator.setInputValue("autoOptimize", true);
        deltaExchangeOperator.initAndValidate();
        deltaExchangeOperator.setID(parameter.getID() + ".deltaExchange");
        adaptableOperatorSampler.setInputValue("operator", deltaExchangeOperator);
        weight += deltaExchangeOperator.getWeight();


        //add Extra DeltaExchangeOperator and SwapOperator
        // heterozygous pairs
        String[] filterArray = {
                "0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0",
                "0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0",
                "0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0",
                "0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0",
                "0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0",
                "0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0"
        };

        IntegerParameter[] pairs0 = new IntegerParameter[filterArray.length];
        BooleanParameter[] pairs1 = new BooleanParameter[filterArray.length];
        for (int i = 0; i < filterArray.length; i++) {
            String filter = filterArray[i];
            pairs0[i] = new IntegerParameter(filter);
            pairs1[i] = new BooleanParameter(filter);
        };


        for (int i = 0; i < pairs0.length; i++) {
            DeltaExchangeOperator extraDeltaExchangeOperator = new DeltaExchangeOperator();
            extraDeltaExchangeOperator.setInputValue("parameter", parameter);
            extraDeltaExchangeOperator.setInputValue("autoOptimize", true);
            extraDeltaExchangeOperator.setInputValue("delta", 1.0 / 16);
            extraDeltaExchangeOperator.setInputValue("weightvector", pairs0[i]);
            extraDeltaExchangeOperator.setInputValue("weight", "2.0"); // set operator weight to 1
            extraDeltaExchangeOperator.setID("deltaExchangePair" + (i + 1));
            extraDeltaExchangeOperator.initAndValidate();
            adaptableOperatorSampler.setInputValue("operator", extraDeltaExchangeOperator);
            weight += extraDeltaExchangeOperator.getWeight();

            SwapOperator swapOperator = new SwapOperator();
            swapOperator.setInputValue("parameter", parameter);
            swapOperator.setInputValue("filter", pairs1[i]);
            swapOperator.setInputValue("weight", "1.0"); // set operator weight to 0.5
            swapOperator.setID("swapPair" + (i + 1));
            swapOperator.initAndValidate();
            adaptableOperatorSampler.setInputValue("operator", swapOperator);
            weight += swapOperator.getWeight();
        }

        //add AdaptableVarianceMultivariateNormalOperator
        AdaptableVarianceMultivariateNormalOperator adaptableVarianceMultivariateNormalOperator = new AdaptableVarianceMultivariateNormalOperator();
        List<Transform> transforms = new ArrayList<>();
        List<Function> logFunctions = new ArrayList<>();
        logFunctions.add(parameter);
        Transform.LogConstrainedSumTransform logConstrainedSumTransform = new Transform.LogConstrainedSumTransform();
        logConstrainedSumTransform.setInputValue("f", logFunctions);
        logConstrainedSumTransform.initAndValidate();
        transforms.add(logConstrainedSumTransform);
        adaptableVarianceMultivariateNormalOperator.setInputValue("weight", weight);
        adaptableVarianceMultivariateNormalOperator.setInputValue("coefficient", 1.0);
        adaptableVarianceMultivariateNormalOperator.setInputValue("scaleFactor", 1.0);
        adaptableVarianceMultivariateNormalOperator.setInputValue("beta", 0.05);
        adaptableVarianceMultivariateNormalOperator.setInputValue("initial", 200 * transforms.size());
        adaptableVarianceMultivariateNormalOperator.setInputValue("burnin", 100 * transforms.size());
        adaptableVarianceMultivariateNormalOperator.setInputValue("every", 1);
        adaptableVarianceMultivariateNormalOperator.setInputValue("allowNonsense", false);
        adaptableVarianceMultivariateNormalOperator.setInputValue("transformations", transforms);
        adaptableVarianceMultivariateNormalOperator.initAndValidate();
        adaptableVarianceMultivariateNormalOperator.setID(parameter.getID() + ".AVMN");
        adaptableOperatorSampler.setInputValue("operator", adaptableVarianceMultivariateNormalOperator);

        adaptableOperatorSampler.setInputValue("weight", weight);
        adaptableOperatorSampler.initAndValidate();
        context.addExtraOperator(adaptableOperatorSampler);
    }

    @Override
    public Class<GT10> getGeneratorClass() { return GT10.class; }

    @Override
    public Class<phylonco.beast.evolution.substitutionmodel.GT10> getBEASTClass() {
        return phylonco.beast.evolution.substitutionmodel.GT10.class;
    }
}
