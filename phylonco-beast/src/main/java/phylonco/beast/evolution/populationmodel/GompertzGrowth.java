package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.operator.UpDownOperator;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("Coalescent intervals for a gompertz growing population.")
public class GompertzGrowth extends PopulationFunction.Abstract implements Loggable, PopFuncWithUpDownOp {
    final public Input<Function> f0Input = new Input<>("f0",
            "Initial proportion of the carrying capacity.", Input.Validate.REQUIRED);
    final public Input<Function> bInput = new Input<>("b",
            "Initial growth rate of tumor growth. Should be greater than 0.", Input.Validate.REQUIRED);
    final public Input<Function> N0Input = new Input<>("N0",
            "Initial population size.", Input.Validate.REQUIRED);

    private AdaptableVarianceMultivariateNormalOperator avmnOperator;

    public GompertzGrowth() {
        // Example of setting up inputs with default values
        // f0Input.setValue(new RealParameter("0.5"), this);
        // bInput.setValue(new RealParameter("0.1"), this);
        // N0Input.setValue(new RealParameter("100"), this);
    }

    @Override
    public void initAndValidate() {
        if (f0Input.get() != null && f0Input.get() instanceof RealParameter) {
            RealParameter f0Param = (RealParameter) f0Input.get();
            f0Param.setBounds(Math.max(0.0, f0Param.getLower()), f0Param.getUpper());
        }

        if (bInput.get() != null && bInput.get() instanceof RealParameter) {
            RealParameter bParam = (RealParameter) bInput.get();
            bParam.setBounds(0.0, Double.POSITIVE_INFINITY);  // b should be positive for Gompertz growth
        }

        if (N0Input.get() != null && N0Input.get() instanceof RealParameter) {
            RealParameter N0Param = (RealParameter) N0Input.get();
            N0Param.setBounds(Math.max(0.0, N0Param.getLower()), N0Param.getUpper());
        }

        // Setup AVMN operator with transformations
        setupAVMNOperator();
    }

    private void setupAVMNOperator() {
        avmnOperator = new AdaptableVarianceMultivariateNormalOperator();
        avmnOperator.setID("AVMNOperator");
        avmnOperator.setInputValue("beta", 0.05);
        avmnOperator.setInputValue("burnin", 400);
        avmnOperator.setInputValue("initial", 800);
        avmnOperator.setInputValue("weight", 2.0);

        List<Transform> transformations = new ArrayList<>();

        // f0 transformation using Interval
        Transform.Interval f0Transform = new Transform.Interval();
        f0Transform.setInputValue("lower", 0.0);
        f0Transform.setInputValue("upper", 1.0);
        f0Transform.setParameter((RealParameter) f0Input.get());
        transformations.add(f0Transform);

        // b transformation using LogTransform
        Transform.LogTransform bTransform = new Transform.LogTransform();
        bTransform.setParameter((RealParameter) bInput.get());
        transformations.add(bTransform);

        avmnOperator.setInputValue("transformations", transformations);
        avmnOperator.initAndValidate();
    }

    public double getNInfinity() {
        return getN0() / getF0();
    }

    public double getF0() {
        return f0Input.get().getArrayValue();
    }

    public final double getGrowthRateB() {
        return bInput.get().getArrayValue();
    }

    public double getN0() {
        return N0Input.get().getArrayValue();
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (f0Input.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) f0Input.get()).getID());
        if (bInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) bInput.get()).getID());
        if (N0Input.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface) N0Input.get()).getID());
        return ids;
    }

    @Override
    public double getPopSize(double t) {
        double f0 = getF0();
        double b = getGrowthRateB();
        double NInfinity = getNInfinity();
        double N0 = getN0();

        return N0 * Math.exp(Math.log(NInfinity / N0) * (1 - Math.exp(b * t)));
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
    }

    @Override
    public void log(long step, PrintStream printStream) {
    }

    @Override
    public void close(PrintStream printStream) {
    }

    @Override
    public UpDownOperator getUpDownOperator(Tree tree) {
        UpDownOperator upDownOperator = new UpDownOperator();
        String idStr = getID() + "Up" + tree.getID() + "DownOperator";
        upDownOperator.setID(idStr);
        upDownOperator.setInputValue("scaleFactor", 0.75);
        upDownOperator.setInputValue("weight", 3.0);
        upDownOperator.setInputValue("up", f0Input.get());
        upDownOperator.setInputValue("up", bInput.get());
        upDownOperator.setInputValue("down", tree);
        upDownOperator.initAndValidate();
        return upDownOperator;
    }
}
