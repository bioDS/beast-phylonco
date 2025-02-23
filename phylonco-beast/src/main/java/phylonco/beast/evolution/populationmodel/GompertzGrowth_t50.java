package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.operator.UpDownOperator;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Coalescent intervals for a Gompertz-growth population (t50-parameterized),
 * with an optional ancestral population size (NA) controlled by I_na (0 or 1).
 *
 * <p>
 * If I_na=0 or NA <= 0, we ignore NA in the Gompertz formula.
 * If I_na=1 and NA>0, the model includes NA as a baseline.
 * The time t50 marks when population size is half of the carrying capacity NInfinity.
 * </p>
 *
 * <p>
 * Unlike the original version, we do not store a fixed "useAncestralPopulation" or "N0"
 * during init. Instead, we recompute logic (like N0) each time in {@code getPopSize(...)}
 * or other methods, similar to GompertzGrowth_f0. This allows I_na to be truly sampled
 * by MCMC, switching formula at runtime.
 * </p>
 */
@Description("Coalescent intervals for a Gompertz (t50) model with an optional NA, dynamically using I_na in getPopSize.")
public class GompertzGrowth_t50 extends PopulationFunction.Abstract implements Loggable,PopFuncWithUpOp, PopFuncWithAVMNOp {

    // -----------------------------------------------------------------------
    // Inputs
    // -----------------------------------------------------------------------
    final public Input<Function> t50Input = new Input<>(
            "t50",
            "Time at which population is half of carrying capacity (NInfinity).",
            Input.Validate.REQUIRED
    );

    final public Input<Function> bInput = new Input<>(
            "b",
            "Gompertz growth rate (> 0).",
            Input.Validate.REQUIRED
    );

    final public Input<Function> NInfinityInput = new Input<>(
            "NInfinity",
            "Carrying capacity of the Gompertz model (must be >0).",
            Input.Validate.REQUIRED
    );

    final public Input<Function> NAInput = new Input<>(
            "NA",
            "Ancestral population size (>=0). If I_na=1 and NA>0, use a baseline offset.",
            Input.Validate.OPTIONAL
    );

    final public Input<IntegerParameter> indicatorParameterInput = new Input<>(
            "I_na",
            "Indicator parameter: 0 or 1. If 1 and NA>0 => use NA; otherwise ignore.",
            Input.Validate.OPTIONAL
    );

    // No local booleans or N0 fields stored: we do "lazy" checking in getPopSize().
    // This ensures MCMC changes to I_na or NA are recognized on each call.

    @Override
    public void initAndValidate() {
        // 1) Basic bounds
        if (t50Input.get() instanceof RealParameter) {
            RealParameter t50Param = (RealParameter) t50Input.get();
            t50Param.setBounds(Math.max(0.0, t50Param.getLower()), t50Param.getUpper());
        }
        if (bInput.get() instanceof RealParameter) {
            RealParameter bParam = (RealParameter) bInput.get();
            bParam.setBounds(Math.max(0.0, bParam.getLower()), bParam.getUpper());
        }
        if (NInfinityInput.get() instanceof RealParameter) {
            RealParameter niParam = (RealParameter) NInfinityInput.get();
            niParam.setBounds(Math.max(0.0, niParam.getLower()), niParam.getUpper());
        }

        // 2) If NA is provided, ensure >=0.
        if (NAInput.get() instanceof RealParameter) {
            RealParameter naParam = (RealParameter) NAInput.get();
            naParam.setBounds(Math.max(0.0, naParam.getLower()), naParam.getUpper());
        }

        // 3) If I_na is provided, clamp it to [0,1].
        if (indicatorParameterInput.get() != null) {
            IntegerParameter iNaParam = indicatorParameterInput.get();
            iNaParam.setInputValue("lower", 0);
            iNaParam.setInputValue("upper", 1);
            int val = iNaParam.getValue();
            if (val < 0) {
                iNaParam.setValue(0);
            } else if (val > 1) {
                iNaParam.setValue(1);
            }
        }
    }

    /**
     * Returns the time t50 (>= 0).
     */
    public double getT50() {
        return t50Input.get().getArrayValue();
    }

    /**
     * Returns the Gompertz growth rate b (>= 0).
     */
    public double getGrowthRate() {
        return bInput.get().getArrayValue();
    }

    /**
     * Returns NInfinity, the carrying capacity.
     */
    public double getNInfinity() {
        return NInfinityInput.get().getArrayValue();
    }

    /**
     * Returns the raw NA from input or 0 if none provided.
     */
    private double getRawNA() {
        if (NAInput.get() != null) {
            return NAInput.get().getArrayValue();
        }
        return 0.0;
    }

    /**
     * Returns the indicator I_na (0 or 1). If not set, we treat as 0.
     */
    private int getI_na() {
        if (indicatorParameterInput.get() == null) {
            return 0;
        }
        return indicatorParameterInput.get().getValue();
    }

    /**
     * Checks if we are effectively using NA: I_na=1 and NA>0.
     */
    private boolean isUsingNA() {
        return (getI_na() == 1 && getRawNA() > 0.0);
    }

    /**
     * Dynamically computes the initial size N0 at t=0 based on t50, b, NInfinity, and NA if used.
     * <p>
     * If using NA:
     *   N0 = (NInfinity - NA) * exp( -ln(2) / exp(b * t50) ) + NA
     * Else:
     *   N0 = NInfinity * 2^(- exp(-b*t50))
     */
    private double computeN0() {
        double t50 = getT50();
        double b   = getGrowthRate();
        double ni  = getNInfinity();
        double na  = getRawNA();

        if (isUsingNA()) {
            double exponent = -Math.log(2) / Math.exp(b * t50);
            return (ni - na) * Math.exp(exponent) + na;
        } else {
            return ni * Math.pow(2, -Math.exp(-b * t50));
        }
    }

    /**
     * Returns the population size N(t). We recalc N0 each time so that changes in I_na or NA
     * (due to MCMC sampling) are reflected immediately.
     */
    @Override
    public double getPopSize(double t) {
        double b   = getGrowthRate();
        double ni  = getNInfinity();
        double na  = getRawNA();
        double N0  = computeN0(); // dynamic

        if (isUsingNA()) {
            // Extended Gompertz with baseline NA
            double N0_minus_na   = N0 - na;
            double ni_minus_na   = ni - na;
            double exponent      = Math.log(ni_minus_na / N0_minus_na)
                    * (1 - Math.exp(b * t));
            return N0_minus_na * Math.exp(exponent) + na;
        } else {
            double logRatio = Math.log(ni / N0);
            double exponent = (1 - Math.exp(b * t)) * logRatio;
            return N0 * Math.exp(exponent);
        }
    }

    /**
     * Numerically integrates ∫(0..t) 1/N(u) du. If I_na flips or NA changes, we get an updated integrand
     * at each step, consistent with the new parameter values.
     */
    @Override
    public double getIntensity(double t) {
        if (t <= 0.0) {
            return 0.0;
        }
        UnivariateFunction integrand = x -> 1.0 / Math.max(getPopSize(x), 1e-20);
        IterativeLegendreGaussIntegrator integrator =
                new IterativeLegendreGaussIntegrator(5, 1.0e-12, 1.0e-8, 2, 10000);
        double intensity = 0.0;
        try {
            intensity = integrator.integrate(100000, integrand, 0.0, t);
        } catch (TooManyEvaluationsException ex) {
            return intensity;
        }
        return intensity;
    }

    /**
     * Not implemented, or we could do a numeric root solver if needed.
     */
    @Override
    public double getInverseIntensity(double x) {
        return 0.0;
    }

    @Override
    public void init(PrintStream out) {
        out.println("# GompertzGrowth_t50 with dynamic I_na usage (f0-like).");
    }

    @Override
    public void log(long step, PrintStream out) {
        double na  = getRawNA();
        int iNaVal = getI_na();
        out.printf("%d\tNInfinity=%.4f\tb=%.4f\tNA=%.4f\tI_na=%d\n",
                step, getNInfinity(), getGrowthRate(), na, iNaVal);
    }

    @Override
    public void close(PrintStream out) {
        out.println("# End GompertzGrowth_t50 log");
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (t50Input.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) t50Input.get()).getID());
        }
        if (bInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) bInput.get()).getID());
        }
        if (NInfinityInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) NInfinityInput.get()).getID());
        }
        if (NAInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) NAInput.get()).getID());
        }
        if (indicatorParameterInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) indicatorParameterInput.get()).getID());
        }
        return ids;
    }

    @Override
    public UpDownOperator getUpOperator(Tree tree) {
        UpDownOperator upDownOperator = new UpDownOperator();
        String idStr = getID() + "Up" + tree.getID() + "DownOperator";
        upDownOperator.setID(idStr);
        upDownOperator.setInputValue("scaleFactor", 0.75);
        upDownOperator.setInputValue("weight", 3.0);
        upDownOperator.setInputValue("up", t50Input.get());
        upDownOperator.setInputValue("up", bInput.get());
        upDownOperator.setInputValue("up", tree);
        upDownOperator.initAndValidate();
        return upDownOperator;
    }

    @Override
    public AdaptableVarianceMultivariateNormalOperator getAVMNOperator(Tree tree) {

        AdaptableVarianceMultivariateNormalOperator avmnOp = new AdaptableVarianceMultivariateNormalOperator();
        String opID = getID() + "AVMNOperator";
        avmnOp.setID(opID);

        avmnOp.setInputValue("beta", 0.05);
        avmnOp.setInputValue("burnin", 400);
        avmnOp.setInputValue("initial", 800);
        avmnOp.setInputValue("weight", 2.0);

        Transform.NoTransform t50Transform = new Transform.NoTransform();
        t50Transform.setID("t50Transform");
        t50Transform.setInputValue("f", t50Input.get());

        t50Transform.initAndValidate();

        Transform.LogTransform bTransform = new Transform.LogTransform();
        bTransform.setID("bTransform");
        bTransform.setInputValue("f", bInput.get());

        bTransform.initAndValidate();

        List<Transform> transforms = new ArrayList<>();
        transforms.add(t50Transform);
        transforms.add(bTransform);

        avmnOp.setInputValue("transformations", transforms);

        avmnOp.initAndValidate();

        return avmnOp;
    }

}
