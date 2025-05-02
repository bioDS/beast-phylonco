package phylonco.lphybeast.tobeast.values;

import NestedBD.evolution.substitutionmodel.BD;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import phylonco.lphy.evolution.copynumbermodel.CopyNumberBD;

public class CopyNumberBDToBEAST implements ValueToBEAST<CopyNumberBD, BD> {

    public static int DEFAULT_NSTATES = 10;

    @Override
    public BD valueToBEAST(Value<CopyNumberBD> value, BEASTContext context) {
        CopyNumberBD copyNumberBD = value.value();

        // Create the BEAST BD model
        BD bdModel = new BD();

        // Set DEFAULT nstates
        RealParameter nstateParam = new RealParameter(String.valueOf(DEFAULT_NSTATES));

        // Get lambda and mu from the LPhy model
        double lambda = copyNumberBD.getLambda().value();
        double mu = copyNumberBD.getMu().value();
        // We assume birth rate (lambda) = death rate (mu)
        double bdRate = lambda;
        // Set the birth-death rate (bdRate) parameter
        RealParameter bdRateParam = new RealParameter(String.valueOf(bdRate));


        // Set inputs of BD model
        bdModel.setInputValue("nstate", nstateParam);
        bdModel.setInputValue("bdRate", bdRateParam);
        bdModel.initAndValidate();

        // Return the completed BD model
        return bdModel;
    }

    @Override
    public Class getValueClass() {
        return CopyNumberBD.class;
    }

    @Override
    public Class<BD> getBEASTClass() {
        return BD.class;
    }
}
