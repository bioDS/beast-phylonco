package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Generator;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.beast.evolution.populationmodel.GompertzGrowth;

public class GompertzToBEAST implements GeneratorToBEAST<GompertzGrowth, phylonco.beast.evolution.populationmodel.GompertzGrowth> {
//     <populationModel id="gompertzPopulationModel" spec="phylonco.beast.evolution.populationmodel.GompertzGrowth">
//                        <parameter name="f0" idref="f0"/>
//                        <parameter name="NInfinity" idref="NInfinity"/>
//                        <parameter name="b" idref="b"/>
//                    </populationModel>

    @Override
    public phylonco.beast.evolution.populationmodel.GompertzGrowth generatorToBEAST(GompertzGrowth gompertzGrowth, BEASTInterface value, BEASTContext context) {
        phylonco.beast.evolution.populationmodel.GompertzGrowth beastGompertzGrowth = new phylonco.beast.evolution.populationmodel.GompertzGrowth();

        double f0 = gompertzGrowth.getF0();


        double b = gompertzGrowth.getGrowthRateB();
        double NInfinity = gompertzGrowth.getNInfinity();

        RealParameter f0Param = new RealParameter(Double.toString(f0));
        RealParameter bParam = new RealParameter(Double.toString(b));
        RealParameter NInfinityParam = new RealParameter(Double.toString(NInfinity));

        beastGompertzGrowth.setInputValue("f0", f0Param);
        beastGompertzGrowth.setInputValue("b", bParam);
        beastGompertzGrowth.setInputValue("NInfinity", NInfinityParam);


        beastGompertzGrowth.initAndValidate();

        return beastGompertzGrowth;
    }



    @Override
    public Class<GompertzGrowth> getGeneratorClass() {
        return GompertzGrowth.class;
    }

    @Override
    public Class<phylonco.beast.evolution.populationmodel.GompertzGrowth> getBEASTClass() {
        return phylonco.beast.evolution.populationmodel.GompertzGrowth.class;
    }

}
