package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Population model with piecewise exponential growth.
 */
@Description("Expansion exponential growth population model without N0 as an independent parameter.")
public class ExpansionGrowth extends PopulationFunction.Abstract implements Loggable {
    final public Input<Function> tauInput = new Input<>("tau",
            "Time before which population size is constant.", Input.Validate.REQUIRED);
    final public Input<Function> rInput = new Input<>("r",
            "The exponential growth rate of the population.", Input.Validate.REQUIRED);
    final public Input<Function> NCInput = new Input<>("NC",
            "Current effective population size after time x.", Input.Validate.REQUIRED);
    final public Input<Function> xInput = new Input<>("x",
            "Transition point time at which growth starts.", Input.Validate.REQUIRED);

    private double N0;

    @Override
    public void initAndValidate() {
        double NC = NCInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double tau = tauInput.get().getArrayValue();
        double x;

        // Loop until we get a valid x within the range [0, tau]
        do {
            x = xInput.get().getArrayValue();

            // Logging or printing a warning to indicate retrying (optional)
            if (x < 0 || x > tau) {
                System.err.println("Warning: Transition time x not in range [0, tau]. Retrying...");
            }

        } while (x < 0 || x > tau); // Repeat if x is outside the range [0, tau]

        if (NC <= 0) {
            throw new IllegalArgumentException("Population size NC must be greater than 0.");
        }
        if (r <= 0) {
            throw new IllegalArgumentException("Growth rate r must be greater than 0.");
        }
        if (tau < 0) {
            throw new IllegalArgumentException("Time tau must be non-negative.");
        }

        // Calculate N0 based on the relationship N0 = NC * exp(-r * (tau - x))
        N0 = NC * Math.exp(-r * (tau - x));
    }


    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();

        if (tauInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) tauInput.get()).getID());
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
        double tau = tauInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double NC = NCInput.get().getArrayValue();
        double x = xInput.get().getArrayValue();

        if (t <= x) {
            return NC;
        } else if (t <= tau) {
            return NC * Math.exp(-r * (t - x));
        } else {
            return N0;
        }
    }

    @Override
    public double getIntensity(double t) {
        double tau = tauInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double NC = NCInput.get().getArrayValue();
        double x = xInput.get().getArrayValue();

        if (t <= x) {
            return t / NC;
        } else if (t <= tau) {
            double firstIntegral = x / NC;
            double secondIntegral = (Math.exp(r * (t - x)) - 1) / (r * NC);
            return firstIntegral + secondIntegral;
        } else {
            double firstIntegral = x / NC;
            double secondIntegral = (Math.exp(r * (tau - x)) - 1) / (r * NC);
            double thirdIntegral = (t - tau) / N0;
            return firstIntegral + secondIntegral + thirdIntegral;
        }
    }

    @Override
    public double getInverseIntensity(double x) {
        // This method implementation can be provided based on model requirements.
        return 0.0;
    }

    @Override
    public void init(PrintStream printStream) {
        printStream.println("# Step\tN0\tComputed_N0\ttau\tr\tNC");
    }

    @Override
    public void log(long step, PrintStream printStream) {
        double NC = NCInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double tau = tauInput.get().getArrayValue();
        double x = xInput.get().getArrayValue();

        // Calculate N0 each time to ensure consistency in logging.
        double computedN0 = NC * Math.exp(-r * (tau - x));

        printStream.println(step + "\t" + computedN0 + "\t" + tau + "\t" + r + "\t" + NC);
    }

    @Override
    public void close(PrintStream printStream) {
        printStream.println("# End of log");
    }

}
