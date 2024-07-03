package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("Coalescent intervals for a constant population.")
public class ConstantGrowth extends PopulationFunction.Abstract implements Loggable {
    final public Input<Function> popSizeParameterInput = new Input<>("N0",
            "present-day population size (defaults to 1.0). ");

    @Override
    public void initAndValidate() {
        if (popSizeParameterInput.get() != null && popSizeParameterInput.get() instanceof RealParameter) {
            RealParameter param = (RealParameter) popSizeParameterInput.get();
            param.setBounds(Math.max(0.0, param.getLower()), param.getUpper());
        }
    }

    /**
     * @return constant population size.
     */
    public double getN0() {
        return popSizeParameterInput.get().getArrayValue();
    }

    /**
     * sets constant population size.
     *
     *
     */
    //    public void setN0(double N0) {
    //        this.N0 = N0;
    //    }

    // Implementation of abstract methods
    @Override
    public double getPopSize(double t) {
        return getN0();
    }

    /**
     * Calculates the integral 1/N(x) dx between start and finish.
     */
    @Override
    public double getIntegral(double start, double finish) {
        return (finish - start) / getN0();
    }

    @Override
    public double getIntensity(double t) {
        return t / getN0();
    }

    @Override
    public double getInverseIntensity(double v) {
        return v * getN0();
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (popSizeParameterInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) popSizeParameterInput.get()).getID());
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
