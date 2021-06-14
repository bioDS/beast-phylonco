package lphybeast.phylonco.generators;

import beast.core.BEASTInterface;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import lphy.evolution.substitutionmodel.GT16;
import lphy.graphicalModel.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

public class GT16ToBEAST implements GeneratorToBEAST<GT16, beast.evolution.substitutionmodel.GT16> {

    @Override
    public beast.evolution.substitutionmodel.GT16 generatorToBEAST(GT16 gt16, BEASTInterface value, BEASTContext context) {

        beast.evolution.substitutionmodel.GT16 beastGT16 = new beast.evolution.substitutionmodel.GT16();

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

        return beastGT16;
    }

    @Override
    public Class<GT16> getGeneratorClass() { return GT16.class; }

    @Override
    public Class<beast.evolution.substitutionmodel.GT16> getBEASTClass() {
        return beast.evolution.substitutionmodel.GT16.class;
    }
}
