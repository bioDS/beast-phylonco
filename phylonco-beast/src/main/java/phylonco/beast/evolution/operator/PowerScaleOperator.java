package phylonco.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Scalable;
import beast.base.inference.StateNode;
import beast.base.inference.operator.kernel.KernelOperator;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Scales several parameters together by powers of one factor: parameter k is multiplied by
 * c^power[k], with log c drawn from a (Bactrian) kernel. The log Hastings ratio is
 * sum_k power[k] * dim[k] * log c.
 *
 * <p>Written for the read-count coverage model, where a cell's coverage has mean alpha t s and
 * variance alpha t s + alpha^2 v s^2: t * c, v * c^2 and every s / c leave every coverage
 * distribution, and so the likelihood, unchanged. The posterior then runs along that ridge (log t,
 * log v and mean log s correlated at +-0.99 on CRC09), where proposals that change one parameter, or
 * scale v by c instead of c^2 (UpDownOperator), are rejected. With parameter="@t @v @s" and
 * power="1 2 -1" every proposal stays on the ridge and only the priors decide acceptance.</p>
 */
@Description("Scales parameters together by powers of one factor (parameter k by c^power[k])")
public class PowerScaleOperator extends KernelOperator {

    final public Input<List<Scalable>> parameterInput = new Input<>("parameter",
            "parameters to scale, each by the factor to its power", new ArrayList<>(), Validate.REQUIRED);
    final public Input<String> powerInput = new Input<>("power",
            "one power per parameter, whitespace separated (e.g. '1 2 -1')", Validate.REQUIRED);
    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "magnitude of the scale factor's kernel", 0.1);
    final public Input<Boolean> optimiseInput = new Input<>("optimise",
            "tune the scale factor for a good acceptance rate (default true)", true);
    final public Input<Double> upperInput = new Input<>("upper", "upper limit of the scale factor's magnitude", 10.0);
    final public Input<Double> lowerInput = new Input<>("lower", "lower limit of the scale factor's magnitude", 0.0);

    private double[] powers;
    private double scaleFactor, upper, lower;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        String[] f = powerInput.get().trim().split("\\s+");
        List<Scalable> params = parameterInput.get();
        if (f.length != params.size()) {
            throw new IllegalArgumentException("PowerScaleOperator: " + params.size() + " parameters but "
                    + f.length + " powers ('" + powerInput.get() + "')");
        }
        powers = new double[f.length];
        for (int k = 0; k < f.length; k++) powers[k] = Double.parseDouble(f[k]);
        scaleFactor = scaleFactorInput.get();
        upper = upperInput.get();
        lower = lowerInput.get();
    }

    @Override
    public double proposal() {
        double scale = kernelDistribution.getScaler(0, Double.NaN, scaleFactor);
        double logHR = 0.0;
        try {
            List<Scalable> params = parameterInput.get();
            for (int k = 0; k < params.size(); k++) {
                Scalable p = (Scalable) ((StateNode) params.get(k)).getCurrentEditable(this);
                // scale(x) multiplies every element by x and returns dim * log x
                logHR += p.scale(Math.pow(scale, powers[k]));
            }
            for (Scalable p : params) {
                if (outsideBounds(p)) return Double.NEGATIVE_INFINITY;
            }
        } catch (Exception e) {
            return Double.NEGATIVE_INFINITY;   // an invalid state, e.g. a bound violated inside scale()
        }
        return logHR;
    }

    private static boolean outsideBounds(Scalable node) {
        if (node instanceof RealScalarParam<?> p) return !p.withinBounds(p.get());
        if (node instanceof RealVectorParam<?> p) {
            for (int i = 0; i < p.size(); i++) if (!p.withinBounds(p.get(i))) return true;
        }
        return false;
    }

    @Override
    public void optimize(double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha) + Math.log(scaleFactor);
            setCoercableParameterValue(Math.exp(delta));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double ratio = Math.max(0.5, Math.min(2.0, prob / getTargetAcceptanceProbability()));
        if (prob < 0.10 || prob > 0.40) {
            return "Try setting scaleFactor to about " + new DecimalFormat("#.###").format(scaleFactor * ratio);
        }
        return "";
    }
}
