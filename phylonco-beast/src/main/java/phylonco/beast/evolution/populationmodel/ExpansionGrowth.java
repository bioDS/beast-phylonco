package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Population model with piecewise exponential growth.
 */
@Description("Expansion exponential growth population model.")
public class ExpansionGrowth extends PopulationFunction.Abstract implements Loggable {
    final public Input<Function> N0Input = new Input<>("N0",
            "Initial population size before tau.", Input.Validate.REQUIRED);
    final public Input<Function> tauInput = new Input<>("tau",
            "Time before which population size is constant at N0.", Input.Validate.REQUIRED);
    final public Input<Function> rInput = new Input<>("r",
            "The exponential growth rate of the population.", Input.Validate.REQUIRED);
    final public Input<Function> NCInput = new Input<>("NC",
            "Current effective infection number (population size) after time x.", Input.Validate.REQUIRED);


    private double x;

    @Override
    public void initAndValidate() {

        double N0 = N0Input.get().getArrayValue();
        double tau = tauInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double NC = NCInput.get().getArrayValue();

        if (N0 <= 0) {
            throw new IllegalArgumentException("Initial population size N0 must be greater than 0.");
        }
        if (NC <= N0) {
            throw new IllegalArgumentException("Current population size NC must be greater than N0.");
        }
        if (r <= 0) {
            throw new IllegalArgumentException("Growth rate r must be greater than 0.");
        }
        if (tau < 0) {
            throw new IllegalArgumentException("Time tau must be non-negative.");
        }

        x = tau + (1.0 / r) * Math.log(N0 / NC);
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();

        if (N0Input.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) N0Input.get()).getID());
        }
        if (tauInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) tauInput.get()).getID());
        }
        if (rInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) rInput.get()).getID());
        }
        if (NCInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) NCInput.get()).getID());
        }
        return ids;
    }

    @Override
    public double getPopSize(double t) {
        double N0 = N0Input.get().getArrayValue();
        double tau = tauInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double NC = NCInput.get().getArrayValue();
        double x = this.x;

        if (t <= x) {
            return NC;
        } else if (t <= tau) {
            return NC * Math.exp(-r * (t-x));
        } else {
            return N0;
        }
    }

    @Override
    public double getIntensity(double t) {
        double N0 = N0Input.get().getArrayValue();
        double tau = tauInput.get().getArrayValue();
        double r = rInput.get().getArrayValue();
        double NC = NCInput.get().getArrayValue();
        double x = this.x;

        if (t <= x) {
            // Case 1: t <= x
            return t / NC;
        } else if (t <= tau) {
            // Case 2: x < t <= tau
            double firstIntegral = x / NC;
            double secondIntegral = (Math.exp(r * (t - x)) - 1) / (r * NC);
            return firstIntegral + secondIntegral;
        } else {
            // Case 3: t > tau
            double firstIntegral = x / NC;
            double secondIntegral = (Math.exp(r * (tau - x)) - 1) / (r * NC);
            double thirdIntegral = (t - tau) / N0;
            return firstIntegral + secondIntegral + thirdIntegral;
        }
    }


    @Override
    public double getInverseIntensity(double x) {

        return 0.0;
    }

    @Override
    public void init(PrintStream printStream) {
        printStream.println("# Step\tN0\ttau\tr\tNC");
    }

    @Override
    public void log(long step, PrintStream printStream) {
        printStream.println(step + "\t" + N0Input.get().getArrayValue() + "\t" +
                tauInput.get().getArrayValue() + "\t" +
                rInput.get().getArrayValue() + "\t" +
                NCInput.get().getArrayValue());
    }

    @Override
    public void close(PrintStream printStream) {
        printStream.println("# End of log");
    }
}

