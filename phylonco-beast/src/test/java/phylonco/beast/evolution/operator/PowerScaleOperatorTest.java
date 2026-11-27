package phylonco.beast.evolution.operator;

import beast.base.inference.State;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * PowerScaleOperator samples the right distribution: with t, v, s lognormal and the operator as the
 * only move, the chain stays on the orbit (t0 c, v0 c^2, s0 / c), where log c is exactly normal
 * (the 1/x of each lognormal density cancels the Jacobian). A wrong Hastings ratio shifts the mean.
 */
public class PowerScaleOperatorTest {

    private static final double M1 = 0.2, S1 = 0.3, M2 = 1.0, S2 = 0.5, M3 = -0.1, S3 = 0.7;

    private static double logLognormal(double x, double m, double sd) {
        double z = (Math.log(x) - m) / sd;
        return -Math.log(x) - 0.5 * z * z;
    }

    private static double logTarget(RealScalarParam<?> t, RealScalarParam<?> v, RealVectorParam<?> s) {
        double lp = logLognormal(t.get(), M1, S1) + logLognormal(v.get(), M2, S2);
        for (int i = 0; i < s.size(); i++) lp += logLognormal(s.get(i), M3, S3);
        return lp;
    }

    @Test
    public void testSamplesTheOrbitDistribution() {
        Randomizer.setSeed(11);
        double t0 = 1.3, v0 = 2.0;
        double[] s0 = {0.8, 1.1, 0.9, 1.4, 0.6, 1.2, 1.0, 0.7, 1.3, 0.95};   // 10 cells: a missing Jacobian shifts the mean by ~1 sd
        RealScalarParam<PositiveReal> t = new RealScalarParam<>(t0, PositiveReal.INSTANCE);
        RealScalarParam<PositiveReal> v = new RealScalarParam<>(v0, PositiveReal.INSTANCE);
        RealVectorParam<PositiveReal> s = new RealVectorParam<>(s0.clone(), PositiveReal.INSTANCE);
        t.setID("t");
        v.setID("v");
        s.setID("s");
        State state = new State();
        state.initByName("stateNode", List.of(t, v, s));
        state.initialise();
        PowerScaleOperator op = new PowerScaleOperator();
        op.initByName("parameter", List.of(t, v, s), "power", "1 2 -1", "scaleFactor", 0.4, "optimise", false,
                "weight", 1.0);

        // the exact distribution of u = log c on the orbit
        int n = s0.length;
        double a1 = Math.log(t0), a2 = Math.log(v0);
        double precision = 1 / (S1 * S1) + 4 / (S2 * S2) + n / (S3 * S3);
        double num = (M1 - a1) / (S1 * S1) + 2 * (M2 - a2) / (S2 * S2);
        for (double x : s0) num += (Math.log(x) - M3) / (S3 * S3);
        double mean = num / precision, sd = 1 / Math.sqrt(precision);

        int steps = 400000, burn = 2000;
        double sum = 0, sum2 = 0;
        int kept = 0;
        double lp = logTarget(t, v, s);
        for (int i = 0; i < steps; i++) {
            state.store(i);
            double logHR = op.proposal();
            boolean accept = false;
            double lpNew = Double.NEGATIVE_INFINITY;
            if (logHR > Double.NEGATIVE_INFINITY) {
                lpNew = logTarget(t, v, s);
                accept = Math.log(Randomizer.nextDouble()) < lpNew - lp + logHR;
            }
            if (accept) {
                lp = lpNew;
            } else {
                state.restore();
            }
            state.setEverythingDirty(false);
            if (i >= burn) {
                double u = Math.log(t.get()) - a1;
                sum += u;
                sum2 += u * u;
                kept++;
            }
            // the orbit: v / t^2 and s * t stay fixed
            if (i % 50000 == 0) {
                assertEquals(v0 / (t0 * t0), v.get() / (t.get() * t.get()), 1e-9);
                assertEquals(s0[2] * t0, s.get(2) * t.get(), 1e-9);
            }
        }
        double m = sum / kept, var = sum2 / kept - m * m;
        assertEquals(mean, m, 0.05 * sd, "mean of log c (exact " + mean + ")");
        assertEquals(sd * sd, var, 0.05 * sd * sd, "variance of log c (exact " + sd * sd + ")");
    }

    @Test
    public void testHastingsRatio() {
        RealScalarParam<PositiveReal> t = new RealScalarParam<>(1.0, PositiveReal.INSTANCE);
        RealVectorParam<PositiveReal> s = new RealVectorParam<>(new double[]{1.0, 2.0, 3.0}, PositiveReal.INSTANCE);
        t.setID("t");
        s.setID("s");
        State state = new State();
        state.initByName("stateNode", List.of(t, s));
        state.initialise();
        PowerScaleOperator op = new PowerScaleOperator();
        op.initByName("parameter", List.of(t, s), "power", "2 -1", "scaleFactor", 0.3, "optimise", false, "weight", 1.0);
        state.store(0);
        double logHR = op.proposal();
        double logC = Math.log(t.get()) / 2;              // t was 1, now c^2
        assertEquals((2 * 1 - 1 * 3) * logC, logHR, 1e-12);
        assertEquals(2.0 * Math.exp(-logC), s.get(1), 1e-12);
    }
}
