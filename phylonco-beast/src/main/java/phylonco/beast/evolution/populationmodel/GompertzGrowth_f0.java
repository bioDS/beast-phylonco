package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.operator.UpDownOperator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Coalescent intervals for a Gompertz growing population with an optional ancestral size NA,
 * controlled by an indicator parameter I_na (0 or 1).
 * <p>
 * If I_na=0, the model ignores NA.
 * If I_na=1 and NA>0, the model uses a baseline ancestral population size in the Gompertz formula.
 */
@Description("Coalescent intervals for a Gompertz growing population with an optional ancestral size NA.")
public class GompertzGrowth_f0 extends PopulationFunction.Abstract implements Loggable,PopFuncWithUpDownOp {

    /**
     * Initial proportion of the carrying capacity.
     */
    public final Input<Function> f0Input = new Input<>(
            "f0",
            "Initial proportion of the carrying capacity.",
            Input.Validate.REQUIRED
    );

    /**
     * Growth rate parameter b, should be >= 0.
     */
    public final Input<Function> bInput = new Input<>(
            "b",
            "Growth rate of the population. Should be greater than 0.",
            Input.Validate.REQUIRED
    );

    /**
     * Initial population size N0.
     */
    public final Input<Function> N0Input = new Input<>(
            "N0",
            "Initial population size.",
            Input.Validate.REQUIRED
    );

    /**
     * Optional ancestral population size (NA). If I_na=1 and NA>0, use it in the extended Gompertz formula.
     */
    public final Input<Function> NAInput = new Input<>(
            "NA",
            "Ancestral population size (baseline). Must be >=0 and <= NInfinity.",
            Input.Validate.OPTIONAL
    );

    /**
     * Optional indicator parameter I_na: 0 or 1. If 1 and NA>0, we incorporate NA; otherwise ignore NA.
     */
    public final Input<IntegerParameter> indicatorParameterInput = new Input<>(
            "I_na",
            "Indicator parameter: 0 or 1. If 1 and NA>0, the ancestral population size will be used.",
            Input.Validate.OPTIONAL
    );

    @Override
    public void initAndValidate() {
        // 1) Ensure f0 >= 0
        if (f0Input.get() instanceof RealParameter) {
            RealParameter f0Param = (RealParameter) f0Input.get();
            f0Param.setBounds(Math.max(0.0, f0Param.getLower()), f0Param.getUpper());
        }

        // 2) Ensure b >= 0
        if (bInput.get() instanceof RealParameter) {
            RealParameter bParam = (RealParameter) bInput.get();
            bParam.setBounds(0.0, Double.POSITIVE_INFINITY);
        }

        // 3) Ensure N0 >= 0
        if (N0Input.get() instanceof RealParameter) {
            RealParameter N0Param = (RealParameter) N0Input.get();
            N0Param.setBounds(Math.max(0.0, N0Param.getLower()), N0Param.getUpper());
        }

        // 4) If NA is provided, ensure 0 <= NA <= NInfinity (where NInfinity = N0/f0)
        double N0 = getN0();
        double f0 = getF0();
        double NInfinity = N0 / f0;  // could be +∞ if f0 is extremely small

        if (NAInput.get() instanceof RealParameter) {
            RealParameter NAParam = (RealParameter) NAInput.get();
            NAParam.setBounds(Math.max(0.0, NAParam.getLower()), Math.min(NInfinity, NAParam.getUpper()));
        }

        // 5) If I_na is provided, clamp it to [0,1].
        if (indicatorParameterInput.get() != null) {
            IntegerParameter iNaParam = indicatorParameterInput.get();
            iNaParam.setInputValue("lower", 0);
            iNaParam.setInputValue("upper", 1);
            int currentValue = iNaParam.getValue();
            if (currentValue < 0) {
                iNaParam.setValue(0);
            } else if (currentValue > 1) {
                iNaParam.setValue(1);
            }
        }
    }

    /**
     * Private helper to retrieve I_na (0 or 1). If not provided, defaults to 0.
     */
    private double getIndicatorValue() {
        if (indicatorParameterInput.get() == null) {
            return 0.0;
        }
        return indicatorParameterInput.get().getValue();
    }

    /**
     * Checks whether NA is actually being used (I_na=1 and NA>0).
     */
    private boolean isUsingNA() {
        return (getIndicatorValue() > 0.0 && getNA() > 0.0);
    }

    /**
     * Computes N_infinity = N0 / f0.
     * This is the carrying capacity that the population will approach as t -> infinity in a Gompertz model.
     *
     * @return NInfinity
     */
    public double getNInfinity() {
        return getN0() / getF0();
    }

    /**
     * Returns the initial proportion f0 (must be >0).
     */
    public double getF0() {
        return f0Input.get().getArrayValue();
    }

    /**
     * Returns the growth rate b.
     */
    public double getGrowthRateB() {
        return bInput.get().getArrayValue();
    }

    /**
     * Returns the initial population size N0.
     */
    public double getN0() {
        return N0Input.get().getArrayValue();
    }

    /**
     * Returns the raw NA if provided, else 0.0.
     */
    private double getRawNA() {
        if (NAInput.get() != null) {
            return NAInput.get().getArrayValue();
        }
        return 0.0;
    }

    /**
     * Returns the effective NA.
     * If I_na=0, or NA <= 0, then 0. Otherwise the raw NA.
     */
    public double getNA() {
        double rawNA = getRawNA();
        if (getIndicatorValue() == 1.0 && rawNA > 0.0) {
            return rawNA;
        }
        return 0.0;
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (f0Input.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) f0Input.get()).getID());
        }
        if (bInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) bInput.get()).getID());
        }
        if (N0Input.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) N0Input.get()).getID());
        }
        if (NAInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) NAInput.get()).getID());
        }
        if (indicatorParameterInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) indicatorParameterInput.get()).getID());
        }
        return ids;
    }

    /**
     * Returns the population size at time t according to the Gompertz model.
     * <p>
     * If isUsingNA() => extended formula with baseline NA:
     * <pre>
     *   N(t) = (N0 - NA)*exp( log((N∞ - NA)/(N0 - NA)) * (1 - exp(b*t)) ) + NA
     * </pre>
     * else => original Gompertz formula:
     * <pre>
     *   N(t) = N0 * exp( log(N∞/N0) * (1 - exp(b*t)) )
     * </pre>
     * @param t The time (>=0).
     * @return  The population size at time t.
     */
    @Override
    public double getPopSize(double t) {
        double b = getGrowthRateB();
        double N0 = getN0();
        double f0 = getF0();
        double NInfinity = N0 / f0;
        double NA = getNA();

        if (isUsingNA()) {
            // Extended Gompertz formula with baseline NA
            double numerator   = NInfinity - NA;
            double denominator = N0 - NA;
            double logRatio    = Math.log(numerator / denominator);
            double exponent    = (1 - Math.exp(b * t)) * logRatio;
            return (N0 - NA) * Math.exp(exponent) + NA;
        } else {
            // Original Gompertz formula
            double logRatio = Math.log(NInfinity / N0);
            double exponent = (1 - Math.exp(b * t)) * logRatio;
            return N0 * Math.exp(exponent);
        }
    }

    /**
     * Calculates the coalescent intensity by integrating 1/N(t) from 0 to t (numerical).
     */
    @Override
    public double getIntensity(double t) {
        if (t == 0) return 0;
        UnivariateFunction function = time -> 1 / Math.max(getPopSize(time), 1e-20);
        IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(5, 1.0e-12, 1.0e-8, 2, 10000);
        double intensity = 0;
        try {
            intensity = integrator.integrate(10000, function, 0, t);
        } catch (TooManyEvaluationsException ex) {
            return intensity;
        }
        return intensity;
    }

    @Override
    public double getInverseIntensity(double x) {
        // Not implemented
        return 0.0;
    }

    @Override
    public void init(PrintStream printStream) {
        // Header for logging
        printStream.println("# Step\tf0\tNInfinity\tb\tNA\tI_na");
    }

    @Override
    public void log(long step, PrintStream printStream) {
        printStream.println(step + "\t" + getF0() + "\t" + getNInfinity()
                + "\t" + getGrowthRateB() + "\t" + getRawNA()
                + "\t" + (int) getIndicatorValue());
    }

    @Override
    public void close(PrintStream printStream) {
        printStream.println("# End of log");
    }

    @Override
    public UpDownOperator getUpDownOperator1(Tree tree) {
        UpDownOperator upDownOperator = new UpDownOperator();
        String idStr = getID() + "Up" + tree.getID() + "DownOperator1";
        upDownOperator.setID(idStr);
        upDownOperator.setInputValue("scaleFactor", 0.75);
        upDownOperator.setInputValue("weight", 3.0);
        upDownOperator.setInputValue("up", f0Input.get());
        upDownOperator.setInputValue("down", tree);
        upDownOperator.initAndValidate();
        return upDownOperator;
    }
    @Override
    public UpDownOperator getUpDownOperator2(Tree tree) {
        UpDownOperator upDownOperator = new UpDownOperator();
        String idStr = getID() + "Up" + tree.getID() + "DownOperator2";
        upDownOperator.setID(idStr);
        upDownOperator.setInputValue("scaleFactor", 0.75);
        upDownOperator.setInputValue("weight", 3.0);
        upDownOperator.setInputValue("up", bInput.get());
        upDownOperator.setInputValue("down", tree);
        upDownOperator.initAndValidate();
        return upDownOperator;
    }


}
