package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Population model with piecewise exponential growth.
 */
@Description("Cons_Exp_Cons exponential growth population model without N0 as an independent parameter.")
public class ExpansionGrowth extends PopulationFunction.Abstract implements Loggable {
    final public Input<Function> NAInput = new Input<>("NA",
            "Ancestral population size after growth.", Input.Validate.REQUIRED);
    final public Input<Function> rInput = new Input<>("r",
            "The exponential growth rate of the population.", Input.Validate.REQUIRED);
    final public Input<Function> NCInput = new Input<>("NC",
            "Current effective population size after time x.", Input.Validate.REQUIRED);
    final public Input<Function> xInput = new Input<>("x",
            "Transition point time at which growth starts.", Input.Validate.REQUIRED);



    @Override
    public void initAndValidate() {
        double NC = NCInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double NA = NAInput.get().getArrayValue();

        // Validate parameter values
        if (NC <= 0) {
            throw new IllegalArgumentException("Population size NC must be greater than 0.");
        }
        if (NA >= NC) {
            throw new IllegalArgumentException("Ancestral population size NA must be less than initial population size NC.");
        }
        if (r <= 0) {
            throw new IllegalArgumentException("Growth rate r must be greater than 0.");
        }
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
        return ids;
    }

    @Override
    public double getPopSize(double t) {
        double NA = NAInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double NC = NCInput.get().getArrayValue();
        double x = xInput.get().getArrayValue();

        if (t < 0) {
            throw new IllegalArgumentException("Time t cannot be negative.");
        }

        if (t <= x) {
            return NC;
        } else {
            return (NC - NA) * Math.exp(-r * (t - x)) + NA;
        }
    }

    @Override
    public double getIntensity(double t) {
        double NA = NAInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double NC = NCInput.get().getArrayValue();
        double x = xInput.get().getArrayValue();

        if (t < 0) return 0.0;

        if (t <= x) {
            return t / NC;
        } else {
            // Integral from 0 to x: ∫(1 / NC) dt = t / NC
            double firstIntegral = x / NC;

            // Integral from x to t: ∫(1 / [(NC - NA)e^{-r(t' - x)} + NA] ) dt
            // No closed-form solution; use numerical integration

            UnivariateFunction integrand = timePoint -> 1.0 / ((NC - NA) * Math.exp(-r * (timePoint - x)) + NA);
            IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(
                    5, 1.0e-12, 1.0e-8, 2, 10000);
            double secondIntegral = integrator.integrate(Integer.MAX_VALUE, integrand, x, t);

            return firstIntegral + secondIntegral;
        }
    }

    @Override
    public double getInverseIntensity(double x) {
        // This method implementation can be provided based on model requirements.
        return 0.0;
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
