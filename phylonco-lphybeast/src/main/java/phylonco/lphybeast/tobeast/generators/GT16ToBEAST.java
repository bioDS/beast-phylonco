package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.inference.operator.DeltaExchangeOperator;
import beast.base.inference.operator.SwapOperator;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.lphy.evolution.substitutionmodel.GT16;

/**
 * This has to create TreeLikelihood.
 * A ~ PhyloCTMC();
 * D ~ ErrorModel(A);
 */
public class GT16ToBEAST implements GeneratorToBEAST<GT16, phylonco.beast.evolution.substitutionmodel.GT16> {

    @Override
    public phylonco.beast.evolution.substitutionmodel.GT16 generatorToBEAST(GT16 gt16, BEASTInterface value, BEASTContext context) {

        phylonco.beast.evolution.substitutionmodel.GT16 beastGT16 = new phylonco.beast.evolution.substitutionmodel.GT16();

        Value<Double[]> rates = gt16.getRates();
        Value<Double[]> freqs = gt16.getFreq();

        RealParameter ratesParameter = (RealParameter)context.getBEASTObject(rates);
        ratesParameter.setInputValue("keys", "AC AG AT CG CT GT");
        ratesParameter.initAndValidate();

        Frequencies freqsParameter = BEASTContext.createBEASTFrequencies(
                (RealParameter) context.getBEASTObject(freqs),
                "0 1 2 3 4 5 6 7 8 9 a b c d e f");
        freqsParameter.initAndValidate();

        beastGT16.setInputValue("nucRates", ratesParameter);
        beastGT16.setInputValue("frequencies", freqsParameter);
        beastGT16.initAndValidate();

        // set operators on frequency and rates
        addDeltaExchangeOperator(freqsParameter.frequenciesInput.get(), context);
        addDeltaExchangeOperator(ratesParameter, context);

        // add extra operators
        addExtraDeltaExchangeOperators(freqsParameter, context);
        addExtraSwapOperators(freqsParameter, context);

        return beastGT16;
    }

    private void addDeltaExchangeOperator(RealParameter parameter, BEASTContext context) {
        DeltaExchangeOperator operator = new DeltaExchangeOperator();
        operator.setInputValue("parameter", parameter);
        operator.setInputValue("weight", context.getOperatorWeight(parameter.getDimension() - 1));
        operator.setInputValue("autoOptimize", false);
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
            operator.setInputValue("autoOptimize", false);
            operator.setInputValue("delta", 1.0 / 16);
            operator.setInputValue("weightvector", pairs[i]);
            operator.setInputValue("weight", "1.0"); // set operator weight to 1
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

    @Override
    public Class<GT16> getGeneratorClass() { return GT16.class; }

    @Override
    public Class<phylonco.beast.evolution.substitutionmodel.GT16> getBEASTClass() {
        return phylonco.beast.evolution.substitutionmodel.GT16.class;
    }
}
