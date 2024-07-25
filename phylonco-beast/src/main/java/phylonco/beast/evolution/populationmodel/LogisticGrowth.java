package phylonco.beast.evolution.populationmodel;


import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("Coalescent intervals for a logistic growing population.")
public class LogisticGrowth extends PopulationFunction.Abstract implements Loggable {
    final public Input<Function> t50Input = new Input<>("t50", "The time at which the population reaches half of its carrying capacity.", Input.Validate.REQUIRED);
    final public Input<Function> nCarryingCapacityInput = new Input<>("nCarryingCapacity", "The carrying capacity of the population.", Input.Validate.REQUIRED);
    final public Input<Function> bInput = new Input<>("b", "The growth rate of the population.", Input.Validate.REQUIRED);

    public LogisticGrowth() {
        // Example of setting up inputs with default values
    }


    @Override
    public void initAndValidate() {
        if (t50Input.get() != null && t50Input.get() instanceof RealParameter) {
            RealParameter t50Param = (RealParameter) t50Input.get();
            t50Param.setBounds(Math.max(0.0, t50Param.getLower()), t50Param.getUpper());
        }

        if (nCarryingCapacityInput.get() != null && nCarryingCapacityInput.get() instanceof RealParameter) {
            RealParameter nCarryingCapacityInputParam = (RealParameter) nCarryingCapacityInput.get();
            nCarryingCapacityInputParam.setBounds(0.0, Double.POSITIVE_INFINITY);
        }

        if (bInput.get() != null && bInput.get() instanceof RealParameter) {
            RealParameter bParam = (RealParameter) bInput.get();
            bParam.setBounds(Math.max(0.0, bParam.getLower()), bParam.getUpper());
        }
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();

        if (t50Input.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) t50Input.get()).getID());
        }

        if (nCarryingCapacityInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) nCarryingCapacityInput.get()).getID());
        }

        if (bInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) bInput.get()).getID());
        }
        return ids;
    }


    public double getNCarryingCapacity() {
        return nCarryingCapacityInput.get().getArrayValue();
    }

    public final double getGrowthRateB() {
        return bInput.get().getArrayValue();
    }

    public double getT50() {
        return t50Input.get().getArrayValue();
    }


    @Override
    public double getPopSize(double t) {
        double b = getGrowthRateB();
        double nCarryingCapacity = getNCarryingCapacity();
        double t50 = getT50();

        return nCarryingCapacity / (1 + Math.exp(b * (t - t50)));
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
        printStream.println("# Step\tf0\tNInfinity\tb");

    }

    @Override
    public void log(long step, PrintStream printStream) {
        printStream.println(step + "\t" + getT50() + "\t" + getNCarryingCapacity() + "\t" + getGrowthRateB());

    }

    @Override
    public void close(PrintStream printStream) {
        printStream.println("# End of log");

    }
}



