package phylonco.beast.evolution.populationmodel;

import beast.base.core.*;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Coalescent intervals for a logistic growing population with an optional NA.
 * <p>
 * The model is:
 * <ul>
 *   <li>If I_na = 0, ignore NA and use logistic growth from 0 up to carrying capacity (nCarryingCapacity).</li>
 *   <li>If I_na = 1 and NA > 0, population transitions from NA at t << t50 up to nCarryingCapacity, following logistic growth.</li>
 *   <li>If I_na = 1 but NA <= 0, effectively ignore NA (same as I_na=0).</li>
 * </ul>
 */
@Description("Coalescent intervals for a logistic growing population with optional NA and a 0/1 indicator I_na.")
public class LogisticGrowth extends PopulationFunction.Abstract implements Loggable {

    // Required parameters
    public final Input<Function> t50Input = new Input<>(
            "t50",
            "The time at which the population reaches half of its carrying capacity.",
            Input.Validate.REQUIRED
    );
    public final Input<Function> nCarryingCapacityInput = new Input<>(
            "nCarryingCapacity",
            "The carrying capacity of the population.",
            Input.Validate.REQUIRED
    );
    public final Input<Function> bInput = new Input<>(
            "b",
            "The logistic growth rate of the population.",
            Input.Validate.REQUIRED
    );

    // Optional parameters
    public final Input<Function> NAInput = new Input<>(
            "NA",
            "The ancestral population size. Must be between 0 and nCarryingCapacity.",
            Input.Validate.OPTIONAL
    );

    /**
     * I_na: 0 or 1. If 0 => ignore NA, If 1 => use NA if NA>0.
     * Optional; defaults to 1 if not provided.
     */
    public final Input<Function> iNaInput = new Input<>(
            "I_na",
            "Indicator for using NA. 0 or 1. If 0 => ignore NA, if 1 => use NA if NA>0.",
            Input.Validate.OPTIONAL
    );

    public LogisticGrowth() {
        // Default constructor
    }

    @Override
    public void initAndValidate() {
        // 1) t50 ≥ 0
        if (t50Input.get() instanceof RealParameter) {
            RealParameter t50Param = (RealParameter) t50Input.get();
            t50Param.setBounds(Math.max(0.0, t50Param.getLower()), t50Param.getUpper());
        }
        // 2) nCarryingCapacity > 0
        if (nCarryingCapacityInput.get() instanceof RealParameter) {
            RealParameter KParam = (RealParameter) nCarryingCapacityInput.get();
            KParam.setBounds(0.0, Double.POSITIVE_INFINITY);
        }
        // 3) b ≥ 0
        if (bInput.get() instanceof RealParameter) {
            RealParameter bParam = (RealParameter) bInput.get();
            bParam.setBounds(Math.max(0.0, bParam.getLower()), bParam.getUpper());
        }
        // 4) NA ∈ [0, nCarryingCapacity]
        if (NAInput.get() instanceof RealParameter) {
            RealParameter NAParam = (RealParameter) NAInput.get();
            // We do not strictly know the upper bound from code alone,
            // but let's attempt to get it from nCarryingCapacity, or fallback to +∞
            double upper = Double.POSITIVE_INFINITY;
            if (nCarryingCapacityInput.get() instanceof RealParameter) {
                RealParameter KParam = (RealParameter) nCarryingCapacityInput.get();
                // getUpper() might be +∞, so this is a best-effort check
                if (!Double.isInfinite(KParam.getUpper())) {
                    upper = KParam.getUpper();
                }
            }
            NAParam.setBounds(0.0, upper);
        }
        // 5) iNa ∈ [0,1] if provided
        if (iNaInput.get() instanceof IntegerParameter) {
            IntegerParameter iNaParam = (IntegerParameter) iNaInput.get();
            iNaParam.setInputValue("lower", 0);
            iNaParam.setInputValue("upper", 1);
            int currentValue = iNaParam.getValue();
            if (currentValue < 0) {
                iNaParam.setValue(0);
            } else if (currentValue > 1) {
                iNaParam.setValue(1);
            }
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
        if (NAInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) NAInput.get()).getID());
        }
        if (iNaInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) iNaInput.get()).getID());
        }
        return ids;
    }

    /** @return logistic carrying capacity K */
    public double getNCarryingCapacity() {
        return nCarryingCapacityInput.get().getArrayValue();
    }

    /** @return logistic growth rate b */
    public double getGrowthRateB() {
        return bInput.get().getArrayValue();
    }

    /** @return logistic midpoint t50 */
    public double getT50() {
        return t50Input.get().getArrayValue();
    }

    /**
     * Returns the raw NA from input if available.
     * Actual usage depends on I_na.
     */
    public double getRawNA() {
        if (NAInput.get() != null) {
            return NAInput.get().getArrayValue();
        }
        return 0.0;
    }

    /**
     * Returns the indicator I_na (0 or 1).
     * If not provided or can't parse, default is 1.
     */
    public double getI_na() {
        if (iNaInput.get() == null) {
            return 1.0;
        }
        return iNaInput.get().getArrayValue();
    }

    /**
     * Returns the effective NA based on I_na logic:
     * <ul>
     *   <li>If I_na=0 => return 0.0</li>
     *   <li>If I_na=1 => return raw NA if >0, else 0.0</li>
     * </ul>
     */
    public double getEffectiveNA() {
        double i_na = getI_na();
        if (i_na == 0.0) {
            // user doesn't want to use NA
            return 0.0;
        }
        double NA = getRawNA();
        return (NA > 0.0) ? NA : 0.0;
    }

    /**
     * The logistic population size function:
     * <pre>
     *   If effective NA > 0:
     *     N(t) = NA + (K - NA) / [1 + exp(b*(t - t50))]
     *   Else
     *     N(t) = K / [1 + exp(b*(t - t50))]
     * </pre>
     */
    @Override
    public double getPopSize(double t) {
        double b = getGrowthRateB();
        double K = getNCarryingCapacity();
        double t50 = getT50();
        double NA = getEffectiveNA();

        if (NA > 0.0) {
            return NA + (K - NA) / (1.0 + Math.exp(b * (t - t50)));
        } else {
            return K / (1.0 + Math.exp(b * (t - t50)));
        }
    }

    /**
     * Returns the coalescent intensity:
     *   Intensity(t) = ∫(0 to t) [1 / N(u)] du
     * via numerical integration.
     */
    @Override
    public double getIntensity(double t) {
        if (t <= 0.0) {
            return 0.0;
        }
        UnivariateFunction integrand = time -> 1.0 / Math.max(getPopSize(time), 1e-20);
        IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(
                5, 1.0e-12, 1.0e-8, 2, 10000
        );
        double intensity = 0.0;
        try {
            intensity = integrator.integrate(100000, integrand, 0, t);
        } catch (TooManyEvaluationsException ex) {
            // If numerical integration fails, return partial result or handle error
            return intensity;
        }
        return intensity;
    }

    /**
     * Inverse intensity is not implemented for logistic growth in this model.
     */
    @Override
    public double getInverseIntensity(double x) {
        return 0.0;
    }

    @Override
    public void init(PrintStream printStream) {
        printStream.println("# step\tt50\tK\tb\tNA\ti_na");
    }

    @Override
    public void log(long step, PrintStream printStream) {
        printStream.println(step + "\t" + getT50()
                + "\t" + getNCarryingCapacity()
                + "\t" + getGrowthRateB()
                + "\t" + getRawNA()
                + "\t" + getI_na());
    }

    @Override
    public void close(PrintStream printStream) {
        printStream.println("# End of logistic growth log");
    }
}
