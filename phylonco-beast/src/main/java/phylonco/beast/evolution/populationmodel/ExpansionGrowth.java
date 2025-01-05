package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.IntegerParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * A piecewise exponential growth model with an optional ancestral population size NA.
 *
 * If I_na=0 or NA<=0, the model ignores NA and behaves like:
 *   N(t) = NC              for 0 <= t <= x
 *   N(t) = NC * exp(-r * (t - x))  for t > x
 *
 * If I_na=1 and NA>0, then for t > x:
 *   N(t) = (NC - NA)*exp(-r*(t - x)) + NA
 *
 * This model has no closed-form intensity when NA>0; a numerical integrator is used for t > x.
 */
@Description("Piecewise exponential growth with optional ancestral size (I_na).")
public class ExpansionGrowth extends PopulationFunction.Abstract implements Loggable {

    // NA: ancestral size used only if I_na=1 and NA>0
    public final Input<Function> NAInput = new Input<>(
            "NA",
            "Ancestral population size (only effective if I_na=1 and NA>0).",
            Input.Validate.REQUIRED
    );

    // r: exponential growth rate
    public final Input<Function> rInput = new Input<>(
            "r",
            "Exponential growth rate (>0).",
            Input.Validate.REQUIRED
    );

    // NC: effective size for [0, x]
    public final Input<Function> NCInput = new Input<>(
            "NC",
            "Population size for t <= x.",
            Input.Validate.REQUIRED
    );

    // x: time boundary where exponential starts
    public final Input<Function> xInput = new Input<>(
            "x",
            "Time at which exponential growth begins.",
            Input.Validate.REQUIRED
    );

    // I_na: indicator (0 or 1). If 1 and NA>0 => use NA, else ignore NA
    public final Input<IntegerParameter> I_naInput = new Input<>(
            "I_na",
            "Indicator for using NA (0 or 1).",
            Input.Validate.OPTIONAL
    );

    @Override
    public void initAndValidate() {
        double NC = NCInput.get().getArrayValue();
        double r  = rInput.get().getArrayValue();
        double NA = NAInput.get().getArrayValue();

        if (NC <= 0.0) {
            throw new IllegalArgumentException("NC must be > 0.");
        }
        if (r <= 0.0) {
            throw new IllegalArgumentException("r must be > 0.");
        }

        // Check I_na validity
        final int iNaVal = getI_naValue();
        if (iNaVal != 0 && iNaVal != 1) {
            throw new IllegalArgumentException("I_na must be 0 or 1.");
        }

        // If I_na=1, then NA>0 must hold, and typically NA < NC in an expansion scenario
        if (iNaVal == 1) {
            if (NA <= 0.0) {
                throw new IllegalArgumentException("NA must be > 0 when I_na=1.");
            }
            if (NA >= NC) {
                throw new IllegalArgumentException("NA must be < NC if I_na=1.");
            }
        }
    }

    /**
     * Returns the integer value of I_na, defaulting to 0 if not provided.
     */
    private int getI_naValue() {
        if (I_naInput.get() == null) {
            return 0; // default to not using NA
        }
        return I_naInput.get().getValue();
    }

    /**
     * Piecewise population size function:
     *   For t <= x: N(t)=NC
     *   For t > x :
     *     if I_na=0 or NA<=0 => N(t)=NC * exp(-r*(t - x))
     *     if I_na=1 and NA>0 => N(t)=(NC - NA)*exp(-r*(t - x)) + NA
     */
    @Override
    public double getPopSize(double t) {
        if (t < 0.0) {
            throw new IllegalArgumentException("Time t cannot be negative.");
        }

        final double NA = NAInput.get().getArrayValue();
        final double r  = rInput.get().getArrayValue();
        final double NC = NCInput.get().getArrayValue();
        final double x  = xInput.get().getArrayValue();
        final int iNa   = getI_naValue();

        if (t <= x) {
            return NC;
        } else {
            if (iNa == 1 && NA > 0.0) {
                // Use ancestral size
                return (NC - NA) * Math.exp(-r * (t - x)) + NA;
            } else {
                // Ignore NA => standard exponential
                return NC * Math.exp(-r * (t - x));
            }
        }
    }

    /**
     * Intensity(t) = ∫ from 0 to t of 1/N(u) du
     * If t <= x => integral = t/NC
     * If t > x :
     *   1) from 0->x => x/NC
     *   2) from x->t => depends on whether NA is used
     *      if I_na=0 => closed form
     *      if I_na=1 => no closed form => numeric integration
     */
    @Override
    public double getIntensity(double t) {
        if (t < 0.0) {
            return 0.0;
        }

        final double NA = NAInput.get().getArrayValue();
        final double r  = rInput.get().getArrayValue();
        final double NC = NCInput.get().getArrayValue();
        final double x  = xInput.get().getArrayValue();
        final int iNa   = getI_naValue();

        if (t <= x) {
            return t / NC;
        } else {
            double firstPart = x / NC; // integral from 0->x

            if (iNa == 1 && NA > 0.0) {
                // Numeric integration from x->t: 1 / [(NC-NA)*exp(-r*(u-x)) + NA]
                UnivariateFunction integrand = timePoint -> 1.0 / ((NC - NA) * Math.exp(-r * (timePoint - x)) + NA);
                IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(
                        5, 1.0e-12, 1.0e-8, 2, 10000);
                double secondPart;
                try {
                    secondPart = integrator.integrate(Integer.MAX_VALUE, integrand, x, t);
                } catch (Exception e) {
                    throw new RuntimeException("Numerical integration failed from x=" + x + " to t=" + t, e);
                }
                return firstPart + secondPart;
            } else {
                // I_na=0 or NA<=0 => closed form: ∫ 1/[NC e^-r(u-x)] = (exp(r*(u-x)))/(r*NC)
                double secondPart = (Math.exp(r * (t - x)) - 1.0) / (r * NC);
                return firstPart + secondPart;
            }
        }
    }

    @Override
    public double getInverseIntensity(double x) {
        // Not implemented
        return 0.0;
    }

    @Override
    public void init(PrintStream printStream) {
        // No initialization needed
    }

    @Override
    public void log(long step, PrintStream printStream) {
        // No logging logic
    }

    @Override
    public void close(PrintStream printStream) {
        // No resources to close
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (NAInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) NAInput.get()).getID());
        }
        if (rInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) rInput.get()).getID());
        }
        if (NCInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) NCInput.get()).getID());
        }
        if (xInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) xInput.get()).getID());
        }
        if (I_naInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) I_naInput.get()).getID());
        }
        return ids;
    }
}
