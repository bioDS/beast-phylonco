package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("Coalescent intervals for a gompertz growing population.")
public class GompertzGrowth_t50 extends PopulationFunction.Abstract implements Loggable {
    final public Input<Function> t50Input = new Input<>("t50",
            "Initial proportion of the carrying capacity.", Input.Validate.REQUIRED);
    final public Input<Function> bInput = new Input<>("b",
            "Initial growth rate of tumor growth. Should be greater than 0.", Input.Validate.REQUIRED);
    final public Input<Function> NInfinityInput = new Input<>("NInfinity",
            "Initial population size.", Input.Validate.REQUIRED);
    final public Input<Function> NAInput = new Input<>("NA",
            "Ancestral population size NA.", Input.Validate.OPTIONAL);

    private AdaptableVarianceMultivariateNormalOperator avmnOperator;



    @Override
    public void initAndValidate() {
        if (t50Input.get() != null && t50Input.get() instanceof RealParameter) {
            RealParameter t50Param = (RealParameter) t50Input.get();
            t50Param.setBounds(Math.max(0.0, t50Param.getLower()), t50Param.getUpper());
        }

        if (bInput.get() != null && bInput.get() instanceof RealParameter) {
            RealParameter bParam = (RealParameter) bInput.get();
            bParam.setBounds(0.0, Double.POSITIVE_INFINITY);
        }

        if (NInfinityInput.get() != null && NInfinityInput.get() instanceof RealParameter) {
            RealParameter NInfinityParam = (RealParameter) NInfinityInput.get();
            NInfinityParam.setBounds(Math.max(0.0, NInfinityParam.getLower()), NInfinityParam.getUpper());
        }
        if (NAInput.get() != null && NAInput.get() instanceof RealParameter) {
            RealParameter NAParam = (RealParameter) NAInput.get();
            NAParam.setBounds(0.0, Double.POSITIVE_INFINITY);  // NA should be non-negative
        }
        if (NAInput.get() != null) {
            double N0 = getN0();
            double NInfinity = getNInfinity();
            double NA = getNA();
            if (N0 <= NA) {
                throw new IllegalArgumentException("Initial population size N0 must be greater than NA.");
            }
            if (NInfinity <= NA) {
                throw new IllegalArgumentException("Carrying capacity NInfinity must be greater than NA.");
            }
        }

    }

//

    public double getNInfinity() {
        return NInfinityInput.get().getArrayValue();
    }


    public final double getGrowthRateB() {
        return bInput.get().getArrayValue();
    }

    public double getT50() {
        return t50Input.get().getArrayValue();
    }
    public double getN0() {
        return getNInfinity() * Math.pow(2, -Math.exp(-getGrowthRateB() * getT50()));
    }

    public double getNA() {
        if (NAInput.get() != null) {
            return NAInput.get().getArrayValue();
        } else {
            return 0.0;
        }
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (t50Input.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) t50Input.get()).getID());
        if (bInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) bInput.get()).getID());
        if (NInfinityInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) NInfinityInput.get()).getID());
        if (NAInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) NAInput.get()).getID());
        return ids;
    }

    @Override
    public double getPopSize(double t) {
        //double t50 = getT50();
        double b = getGrowthRateB();
        double NInfinity = getNInfinity();
        double N0 = getN0();

        return N0 * Math.exp(Math.log(NInfinity / N0) * (1 - Math.exp(b * t)));
    }

    @Override
    public double getIntensity(double t) {
        if (t == 0) return 0;
        UnivariateFunction function = time -> 1 / Math.max(getPopSize(time), 1e-20);
        IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(5, 1.0e-12, 1.0e-8, 2, 10000);
        double intensity = 0;
        try {
            intensity = integrator.integrate(100000, function, 0, t);
        } catch (TooManyEvaluationsException ex) {
            return intensity;
        }
        return intensity;
    }

    @Override
    public double getInverseIntensity(double x) {
        return 0;
    }

    @Override
    public void init(PrintStream printStream) {
    }

    @Override
    public void log(long step, PrintStream printStream) {
    }

    @Override
    public void close(PrintStream printStream) {
    }

//    @Override
//    public UpDownOperator getUpDownOperator1(Tree tree) {
//        UpDownOperator upDownOperator = new UpDownOperator();
//        String idStr = getID() + "Up" + tree.getID() + "DownOperator1";
//        upDownOperator.setID(idStr);
//        upDownOperator.setInputValue("scaleFactor", 0.75);
//        upDownOperator.setInputValue("weight", 3.0);
//        upDownOperator.setInputValue("up", f0Input.get());
//        upDownOperator.setInputValue("down", tree);
//        upDownOperator.initAndValidate();
//        return upDownOperator;
//    }
//
//    @Override
//    public UpDownOperator getUpDownOperator2(Tree tree) {
//        UpDownOperator upDownOperator = new UpDownOperator();
//        String idStr = getID() + "Up" + tree.getID() + "DownOperator2";
//        upDownOperator.setID(idStr);
//        upDownOperator.setInputValue("scaleFactor", 0.75);
//        upDownOperator.setInputValue("weight", 3.0);
//        upDownOperator.setInputValue("up", bInput.get());
//        upDownOperator.setInputValue("down", tree);
//        upDownOperator.initAndValidate();
//        return upDownOperator;
//    }

}
