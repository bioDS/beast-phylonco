package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("Coalescent intervals for a exponentially growing population.")
public class ExponentialGrowth extends PopulationFunction.Abstract implements Loggable {
    final public Input<Function> popSizeParameterInput = new Input<>("N0",
            "present-day population size (defaults to 1.0). ");
    final public Input<Function> growthRateParameterInput = new Input<>("GrowthRate",
            "Growth rate is the exponent of the exponential growth. " +
                    "A value of zero represents a constant population size, negative values represent " +
                    "decline towards the present, positive numbers represents exponential growth towards " +
                    "the present.");


    @Override
    public void initAndValidate() {
        if (popSizeParameterInput.get() != null && popSizeParameterInput.get() instanceof RealParameter) {
            RealParameter param = (RealParameter) popSizeParameterInput.get();
            param.setBounds( Math.max(0.0, param.getLower()), param.getUpper());
        }
//        if (growthRateParameter.get() != null) {
//            growthRateParameter.get().setBounds(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//        }
    }

    /**
     * @return initial population size.
     */
    public double getN0() {
        return popSizeParameterInput.get().getArrayValue();
    }


    /**
     * sets initial population size.
     *
     * @param N0 new size
     */
//    public void setN0(double N0) {
//        this.N0 = N0;
//    }

    /**
     * @return growth rate.
     */
    public final double getGrowthRate() {
        return growthRateParameterInput.get().getArrayValue();
    }

    /**
     * sets growth rate to r.
     *
     * @param r
     */
//    public void setGrowthRate(double r) {
//        this.r = r;
//    }

    /**
     * An alternative parameterization of this model. This
     * function sets growth rate for a given doubling time.
     *
     * @param
     */
//    public void setDoublingTime(double doublingTime) {
//        setGrowthRate(Math.log(2) / doublingTime);
//    }

    // Implementation of abstract methods
    @Override
    public double getPopSize(double t) {

        double r = getGrowthRate();
        if (r == 0) {
            return getN0();
        } else {
            return getN0() * Math.exp(-t * r);
        }
    }

    /**
     * Calculates the integral 1/N(x) dx between start and finish.
     */
    @Override
    public double getIntegral(double start, double finish) {
        double r = getGrowthRate();
        if (r == 0.0) {
            return (finish - start) / getN0();
        } else {
            return (Math.exp(finish * r) - Math.exp(start * r)) / getN0() / r;
        }
    }

    @Override
    public double getIntensity(double t) {
        double r = getGrowthRate();
        if (r == 0.0) {
            return t / getN0();
        } else {
            return (Math.exp(t * r) - 1.0) / getN0() / r;
        }
    }

    @Override
    public double getInverseIntensity(double v) {
        return 0;
    }

//    @Override
//    public double getInverseIntensity(double x) {
//
//        double r = getGrowthRate();
//        if (r == 0.0) {
//            return getN0() * x;
//        } else {
//            return Math.log(1.0 + getN0() * x * r) / r;
//        }
//    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (popSizeParameterInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)popSizeParameterInput.get()).getID());
        return ids;
    }

    @Override
    public void init(PrintStream printStream) {

    }

    @Override
    public void log(long l, PrintStream printStream) {

    }

    @Override
    public void close(PrintStream printStream) {

    }
}




