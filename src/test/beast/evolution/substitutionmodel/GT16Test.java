package test.beast.evolution.substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GT16;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

public class GT16Test {

    private static double DELTA = 1e-10;

    private GT16 gt16;
    private int nrOfStates = 16;

    public void setupModel(Double[] pi, Double[] rates) {
        gt16 = new GT16();

        RealParameter f = new RealParameter(pi);

        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f, "estimate", false);

        gt16.initByName(
                "rateAC", new RealParameter(rates[0].toString()),
                "rateAG", new RealParameter(rates[1].toString()),
                "rateAT", new RealParameter(rates[2].toString()),
                "rateCG", new RealParameter(rates[3].toString()),
                "rateCT", new RealParameter(rates[4].toString()),
                "rateGT", new RealParameter(rates[5].toString()),
                "frequencies", freqs
        );
        nrOfStates = gt16.getStateCount();
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=12)
     * t <- 10
     *
     * rates <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
     * freqs <- rep(1/16, 16)
     *
     * rateAC <- rates[1]
     * rateAG <- rates[2]
     * rateAT <- rates[3]
     * rateCG <- rates[4]
     * rateCT <- rates[5]
     * rateGT <- rates[6]
     *
     * Q <- matrix(c(
     * 0, rateAC, rateAG, rateAT, rateAC, 0, 0, 0, rateAG, 0, 0, 0, rateAT, 0, 0, 0,
     * rateAC, 0, rateCG, rateCT, 0, rateAC, 0, 0, 0, rateAG, 0, 0, 0, rateAT, 0, 0,
     * rateAG, rateCG, 0, rateGT, 0, 0, rateAC, 0, 0, 0, rateAG, 0, 0, 0, rateAT, 0,
     * rateAT, rateCT, rateGT, 0, 0, 0, 0, rateAC, 0, 0, 0, rateAG, 0, 0, 0, rateAT,
     * rateAC, 0, 0, 0, 0, rateAC, rateAG, rateAT, rateCG, 0, 0, 0, rateCT, 0, 0, 0,
     * 0, rateAC, 0, 0, rateAC, 0, rateCG, rateCT, 0, rateCG, 0, 0, 0, rateCT, 0, 0,
     * 0, 0, rateAC, 0, rateAG, rateCG, 0, rateGT, 0, 0, rateCG, 0, 0, 0, rateCT, 0,
     * 0, 0, 0, rateAC, rateAT, rateCT, rateGT, 0, 0, 0, 0, rateCG, 0, 0, 0, rateCT,
     * rateAG, 0, 0, 0, rateCG, 0, 0, 0, 0, rateAC, rateAG, rateAT, rateGT, 0, 0, 0,
     * 0, rateAG, 0, 0, 0, rateCG, 0, 0, rateAC, 0, rateCG, rateCT, 0, rateGT, 0, 0,
     * 0, 0, rateAG, 0, 0, 0, rateCG, 0, rateAG, rateCG, 0, rateGT, 0, 0, rateGT, 0,
     * 0, 0, 0, rateAG, 0, 0, 0, rateCG, rateAT, rateCT, rateGT, 0, 0, 0, 0, rateGT,
     * rateAT, 0, 0, 0, rateCT, 0, 0, 0, rateGT, 0, 0, 0, 0, rateAC, rateAG, rateAT,
     * 0, rateAT, 0, 0, 0, rateCT, 0, 0, 0, rateGT, 0, 0, rateAC, 0, rateCG, rateCT,
     * 0, 0, rateAT, 0, 0, 0, rateCT, 0, 0, 0, rateGT, 0, rateAG, rateCG, 0, rateGT,
     * 0, 0, 0, rateAT, 0, 0, 0, rateCT, 0, 0, 0, rateGT, rateAT, rateCT, rateGT, 0
     * ), nrow=16, byrow=T)
     *
     * d <- -1 * rowSums(Q)
     * diag(Q) <- d
     * beta <- as.vector(-1 / (freqs %*% d))
     * expm(beta * Q * t)
     *
     */
    @Test
    public void testTransitionLong() {
        double t = 10;

        Double p = 1.0 / nrOfStates;

        Double[] pi = {p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p};

        // rates AC, AG, AT, CG, CT, GT
        Double[] rates = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

//        setupModel(pi, rates);

        double[] expected = new double[] {
            0.07211852, 0.06463233, 0.06567749, 0.06612057, 0.06463233, 0.05792324, 0.05885991, 0.05925699, 0.06567749, 0.05885991, 0.05981172, 0.06021523, 0.06612057, 0.05925699, 0.06021523, 0.06062146,
            0.06463233, 0.06853995, 0.06776782, 0.06760882, 0.05792324, 0.06142523, 0.06073325, 0.06059075, 0.05885991, 0.06241853, 0.06171536, 0.06157056, 0.05925699, 0.06283963, 0.06213171, 0.06198593,
            0.06567749, 0.06776782, 0.06765379, 0.06744982, 0.05885991, 0.06073325, 0.06063106, 0.06044826, 0.05981172, 0.06171536, 0.06161151, 0.06142576, 0.06021523, 0.06213171, 0.06202717, 0.06184016,
            0.06612057, 0.06760882, 0.06744982, 0.06736971, 0.05925699, 0.06059075, 0.06044826, 0.06037646, 0.06021523, 0.06157056, 0.06142576, 0.06135280, 0.06062146, 0.06198593, 0.06184016, 0.06176671,
            0.06463233, 0.05792324, 0.05885991, 0.05925699, 0.06853995, 0.06142523, 0.06241853, 0.06283963, 0.06776782, 0.06073325, 0.06171536, 0.06213171, 0.06760882, 0.06059075, 0.06157056, 0.06198593,
            0.05792324, 0.06142523, 0.06073325, 0.06059075, 0.06142523, 0.06513896, 0.06440513, 0.06425402, 0.06073325, 0.06440513, 0.06367958, 0.06353017, 0.06059075, 0.06425402, 0.06353017, 0.06338111,
            0.05885991, 0.06073325, 0.06063106, 0.06044826, 0.06241853, 0.06440513, 0.06429677, 0.06410291, 0.06171536, 0.06367958, 0.06357243, 0.06338076, 0.06157056, 0.06353017, 0.06342327, 0.06323206,
            0.05925699, 0.06059075, 0.06044826, 0.06037646, 0.06283963, 0.06425402, 0.06410291, 0.06402678, 0.06213171, 0.06353017, 0.06338076, 0.06330549, 0.06198593, 0.06338111, 0.06323206, 0.06315696,
            0.06567749, 0.05885991, 0.05981172, 0.06021523, 0.06776782, 0.06073325, 0.06171536, 0.06213171, 0.06765379, 0.06063106, 0.06161151, 0.06202717, 0.06744982, 0.06044826, 0.06142576, 0.06184016,
            0.05885991, 0.06241853, 0.06171536, 0.06157056, 0.06073325, 0.06440513, 0.06367958, 0.06353017, 0.06063106, 0.06429677, 0.06357243, 0.06342327, 0.06044826, 0.06410291, 0.06338076, 0.06323206,
            0.05981172, 0.06171536, 0.06161151, 0.06142576, 0.06171536, 0.06367958, 0.06357243, 0.06338076, 0.06161151, 0.06357243, 0.06346546, 0.06327412, 0.06142576, 0.06338076, 0.06327412, 0.06308335,
            0.06021523, 0.06157056, 0.06142576, 0.06135280, 0.06213171, 0.06353017, 0.06338076, 0.06330549, 0.06202717, 0.06342327, 0.06327412, 0.06319897, 0.06184016, 0.06323206, 0.06308335, 0.06300843,
            0.06612057, 0.05925699, 0.06021523, 0.06062146, 0.06760882, 0.06059075, 0.06157056, 0.06198593, 0.06744982, 0.06044826, 0.06142576, 0.06184016, 0.06736971, 0.06037646, 0.06135280, 0.06176671,
            0.05925699, 0.06283963, 0.06213171, 0.06198593, 0.06059075, 0.06425402, 0.06353017, 0.06338111, 0.06044826, 0.06410291, 0.06338076, 0.06323206, 0.06037646, 0.06402678, 0.06330549, 0.06315696,
            0.06021523, 0.06213171, 0.06202717, 0.06184016, 0.06157056, 0.06353017, 0.06342327, 0.06323206, 0.06142576, 0.06338076, 0.06327412, 0.06308335, 0.06135280, 0.06330549, 0.06319897, 0.06300843,
            0.06062146, 0.06198593, 0.06184016, 0.06176671, 0.06198593, 0.06338111, 0.06323206, 0.06315696, 0.06184016, 0.06323206, 0.06308335, 0.06300843, 0.06176671, 0.06315696, 0.06300843, 0.06293359
        };

        double[] observed = new double[nrOfStates * nrOfStates];

//        model.getTransitionProbabilities(null, t, 0, 1, observed);
//        assertArrayEquals(expected, observed, DELTA);
    }

}
