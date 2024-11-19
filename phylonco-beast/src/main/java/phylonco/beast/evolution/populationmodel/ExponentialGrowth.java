package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

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
    final public Input<Function> ancestralPopulationParameterInput = new Input<>("NA",
            "Ancestral population size. If set, the population size approaches NA as time increases.");



    @Override
    public void initAndValidate() {
        if (popSizeParameterInput.get() != null && popSizeParameterInput.get() instanceof RealParameter) {
            RealParameter param = (RealParameter) popSizeParameterInput.get();
            param.setBounds( Math.max(0.0, param.getLower()), param.getUpper());
        }

        if (ancestralPopulationParameterInput.get() != null && ancestralPopulationParameterInput.get() instanceof RealParameter) {
            RealParameter param = (RealParameter) ancestralPopulationParameterInput.get();
            param.setBounds(Math.max(0.0, param.getLower()), param.getUpper());
        }

    }

    /**
     * @return initial population size.
     */
    public double getN0() {
        return popSizeParameterInput.get().getArrayValue();
    }

    public double getNA() {
        if (ancestralPopulationParameterInput.get() != null) {
            return ancestralPopulationParameterInput.get().getArrayValue();
        } else {
            return 0.0;
        }
    }

    /**
     * @return growth rate.
     */
    public final double getGrowthRate() {
        return growthRateParameterInput.get().getArrayValue();
    }

    public boolean isUsingNA() {
        return ancestralPopulationParameterInput.get() != null && getNA() > 0.0;
    }


    // Implementation of abstract methods
    @Override
    public double getPopSize(double t) {
        double r = getGrowthRate();
        double N0 = getN0();
        double NA = getNA();

        if (r == 0.0) {
            if (isUsingNA()) {
                return NA;
            } else {
                return N0;
            }
        } else {
            if (isUsingNA()) {
                return (N0 - NA) * Math.exp(-r * t) + NA;
            } else {
                return N0 * Math.exp(-r * t);
            }
        }
    }

//    /**
//     * Calculates the integral 1/N(x) dx between start and finish.
//     */
//    @Override
//    public double getIntegral(double start, double finish) {
//        double r = getGrowthRate();
//        if (r == 0.0) {
//            return (finish - start) / getN0();
//        } else {
//            return (Math.exp(finish * r) - Math.exp(start * r)) / getN0() / r;
//        }
//    }

    @Override
    public double getIntensity(double t) {
        if (t == 0.0) {
            return 0.0;
        }

        if (isUsingNA()) {

            UnivariateFunction function = time -> 1 / Math.max(getPopSize(time), 1e-20);
            IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(5, 1.0e-12, 1.0e-8, 2, 10000);

            try {
                return integrator.integrate(Integer.MAX_VALUE, function, 0.0, t);
            } catch (Exception e) {
                throw new RuntimeException("Numerical integration failed for t = " + t, e);
            }
        } else {
            double r = getGrowthRate();
            double N0 = getN0();

            if (r == 0.0) {
                return t / N0;
            } else {
                return (Math.exp(r * t) - 1.0) / (r * N0);
            }
        }
    }

    @Override
    public double getInverseIntensity(double v) {
        return 0;
    }


    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (popSizeParameterInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)popSizeParameterInput.get()).getID());
        if (growthRateParameterInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) growthRateParameterInput.get()).getID());
        if (ancestralPopulationParameterInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) ancestralPopulationParameterInput.get()).getID());
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




