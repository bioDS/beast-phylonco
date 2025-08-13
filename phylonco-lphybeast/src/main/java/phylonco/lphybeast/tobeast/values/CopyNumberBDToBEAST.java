package phylonco.lphybeast.tobeast.values;

import NestedBD.evolution.substitutionmodel.BD;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import phylonco.lphy.evolution.copynumbermodel.CopyNumberBD;

public class CopyNumberBDToBEAST implements ValueToBEAST<CopyNumberBD, BD> {

    public static int DEFAULT_NSTATES = 15;

    @Override
    public BD valueToBEAST(Value<CopyNumberBD> value, BEASTContext context) {
        CopyNumberBD copyNumberBD = value.value();

        // Create the BEAST BD model
        BD bdModel = new BD();

        // Set DEFAULT nstates
        RealParameter nstateParam = new RealParameter(String.valueOf(DEFAULT_NSTATES));

        // Set inputs of BD model
        bdModel.setInputValue("nstate", nstateParam);
//        bdModel.setInputValue("bdRate", bdRateParam);
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
