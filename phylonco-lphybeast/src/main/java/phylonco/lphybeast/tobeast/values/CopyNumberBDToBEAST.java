package phylonco.lphybeast.tobeast.values;

import NestedBD.evolution.substitutionmodel.BD;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import phylonco.lphy.evolution.copynumbermodel.CopyNumberBD;

public class CopyNumberBDToBEAST implements ValueToBEAST<CopyNumberBD, BD> {

    @Override
    public BD valueToBEAST(Value<CopyNumberBD> value, BEASTContext context) {
        CopyNumberBD copyNumberBD = value.value();

        // Create the BEAST BD model
        BD bdModel = new BD();

        // Get nstates
        int nstates = copyNumberBD.getNstate().value();
        RealParameter nstateParam = new RealParameter(String.valueOf(nstates));

        // Get bdRate from the LPhy model
        double bdRate = copyNumberBD.getBdRate().value();
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
