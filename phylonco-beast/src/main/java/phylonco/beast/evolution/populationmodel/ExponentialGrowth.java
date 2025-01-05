package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * A population model representing exponential growth dynamics with an optional ancestral population size.
 * <p>
 * The model can work in two modes depending on the indicator I_na:
 * <ul>
 *   <li>If I_na = 0: The model does not use an ancestral population size (NA). The population size evolves as N(t) = N0 * exp(-r * t).</li>
 *   <li>If I_na = 1 and NA > 0: The model incorporates an ancestral population size NA, resulting in a population size function:
 *       N(t) = (N0 - NA) * exp(-r * t) + NA.</li>
 * </ul>
 * When the growth rate r = 0, the population size is constant (either N0 if I_na=0 or NA if I_na=1 and NA>0).
 *
 * The class uses a numerical integrator (IterativeLegendreGaussIntegrator) for the intensity calculation if NA is used
 * and the model has no closed-form integral solution.
 */
@Description("Coalescent intervals for an exponentially growing population with an optional ancestral population size parameter.")
public class ExponentialGrowth extends PopulationFunction.Abstract implements Loggable {

    // Inputs for parameters
    final public Input<Function> popSizeParameterInput = new Input<>(
            "N0",
            "Present-day population size. Defaults to 1.0 if not specified.",
            Input.Validate.OPTIONAL
    );

    final public Input<Function> growthRateParameterInput = new Input<>(
            "GrowthRate",
            "Growth rate of the exponential model. A value of zero means constant population size.",
            Input.Validate.OPTIONAL
    );

    final public Input<Function> ancestralPopulationParameterInput = new Input<>(
            "NA",
            "Ancestral population size. Used only if I_na=1 and NA>0.",
            Input.Validate.OPTIONAL
    );

    final public Input<IntegerParameter> indicatorParameterInput = new Input<>(
            "I_na",
            "Indicator parameter: 0 or 1. If I_na=1 and NA>0, the ancestral population size will be used.",
            Input.Validate.OPTIONAL
    );

    @Override
    public void initAndValidate() {
        // Ensure N0 is non-negative
        if (popSizeParameterInput.get() instanceof RealParameter) {
            RealParameter n0Param = (RealParameter) popSizeParameterInput.get();
            n0Param.setBounds(Math.max(0.0, n0Param.getLower()), n0Param.getUpper());
        }

        // Ensure NA is non-negative
        if (ancestralPopulationParameterInput.get() instanceof RealParameter) {
            RealParameter naParam = (RealParameter) ancestralPopulationParameterInput.get();
            naParam.setBounds(Math.max(0.0, naParam.getLower()), naParam.getUpper());
        }

        // Validate I_na parameter if provided
        if (indicatorParameterInput.get() != null) {
            IntegerParameter iNaParam = indicatorParameterInput.get();
            final int lowerBound = 0;
            final int upperBound = 1;

            // Set the bounds and clamp current value
            iNaParam.setInputValue("lower", lowerBound);
            iNaParam.setInputValue("upper", upperBound);

            int currentValue = iNaParam.getValue();
            if (currentValue < lowerBound) {
                iNaParam.setValue(lowerBound);
            } else if (currentValue > upperBound) {
                iNaParam.setValue(upperBound);
            }
        }
    }

    /**
     * Retrieves the indicator I_na.
     * @return 0.0 or 1.0 depending on I_na parameter; if not set, returns 0.0.
     */
    private double getIndicatorValue() {
        if (indicatorParameterInput.get() == null) {
            return 0.0;
        }
        return indicatorParameterInput.get().getValue();
    }

    /**
     * Present-day population size N0.
     */
    public double getN0() {
        return popSizeParameterInput.get() != null ? popSizeParameterInput.get().getArrayValue() : 1.0;
    }

    /**
     * Ancestral population size NA.
     * If I_na=1, the model may use NA if NA > 0.
     */
    public double getNA() {
        double naValue = 0.0;
        if (ancestralPopulationParameterInput.get() != null) {
            naValue = ancestralPopulationParameterInput.get().getArrayValue();
        }
        return getIndicatorValue() * naValue;
    }

    /**
     * Exponential growth rate.
     */
    public double getGrowthRate() {
        return growthRateParameterInput.get() != null ? growthRateParameterInput.get().getArrayValue() : 0.0;
    }

    /**
     * Checks if the ancestral population size should be used.
     */
    public boolean isUsingNA() {
        return getNA() > 0.0;
    }

    /**
     * Returns the population size at time t.
     * @param t time (>=0)
     * @return population size at time t
     */
    @Override
    public double getPopSize(double t) {
        final double r = getGrowthRate();
        final double N0 = getN0();
        final double NA = getNA();

        // Handle the constant case
        if (r == 0.0) {
            return isUsingNA() ? NA : N0;
        }

        // When using NA, population size transitions from N0 at t=0 towards NA as t increases.
        if (isUsingNA()) {
            // N(t) = (N0 - NA)*exp(-r*t) + NA
            return (N0 - NA) * Math.exp(-r * t) + NA;
        } else {
            // Standard exponential decay: N(t) = N0*exp(-r*t)
            return N0 * Math.exp(-r * t);
        }
    }

    /**
     * Returns the coalescent intensity at time t:
     * Intensity(t) = ∫(0 to t) (1/N(u)) du.
     * If NA is used and there's no closed form, perform numerical integration.
     * @param t time
     * @return intensity at time t
     */
    @Override
    public double getIntensity(double t) {
        if (t <= 0.0) {
            return 0.0;
        }

        // If using NA, we may need numerical integration due to the piecewise nature and lack of a closed form.
        if (isUsingNA()) {
            UnivariateFunction integrand = time -> 1.0 / Math.max(getPopSize(time), 1e-20);
            IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(
                    5, 1.0e-12, 1.0e-8, 2, 10000
            );
            try {
                return integrator.integrate(Integer.MAX_VALUE, integrand, 0.0, t);
            } catch (Exception e) {
                throw new RuntimeException("Numerical integration failed at t = " + t, e);
            }
        } else {
            // If not using NA, we have simpler formulas:
            final double r = getGrowthRate();
            final double N0 = getN0();

            if (r == 0.0) {
                // Integral of 1/N0 from 0 to t = t/N0
                return t / N0;
            } else {
                // Integral of exp(r*u)/(N0) from 0 to t = (exp(r*t)-1)/(r*N0)
                return (Math.exp(r * t) - 1.0) / (r * N0);
            }
        }
    }

    /**
     * Inverse intensity is not implemented in this model.
     */
    @Override
    public double getInverseIntensity(double v) {
        // Not implemented
        return 0;
    }

    /**
     * Returns a list of parameter IDs involved in this model.
     */
    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (popSizeParameterInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) popSizeParameterInput.get()).getID());
        }
        if (growthRateParameterInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) growthRateParameterInput.get()).getID());
        }
        if (ancestralPopulationParameterInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) ancestralPopulationParameterInput.get()).getID());
        }
        if (indicatorParameterInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) indicatorParameterInput.get()).getID());
        }
        return ids;
    }

    @Override
    public void init(PrintStream printStream) {
        // No initialization required
    }

    @Override
    public void log(long l, PrintStream printStream) {
        // No logging performed by this model
    }

    @Override
    public void close(PrintStream printStream) {
        // No resources to close
    }
}
