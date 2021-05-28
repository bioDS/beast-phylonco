package test.beast.evolution.substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GT16;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

public class GT16Test {

    private static double DELTA = 1e-7;

    private GT16 gt16;
    private int nrOfStates = 16;

    public void setupModel(Double[] pi, Double[] rates) {
        gt16 = new GT16();

        RealParameter f = new RealParameter(pi);

        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f, "estimate", false);

        RealParameter nucRates = new RealParameter(rates);
        nucRates.setInputValue("keys", "AC AG AT CG CT GT");
        nucRates.initAndValidate();
//        for (int i = 0; i < 6; i++) {
//            nucRates.setValue(i, rates[i]);
//        }

        gt16.initByName(
                "nucRates", nucRates,
                "frequencies", freqs
        );
        /*
        gt16.initByName(
                "rateAC", new RealParameter(rates[0].toString()),
                "rateAG", new RealParameter(rates[1].toString()),
                "rateAT", new RealParameter(rates[2].toString()),
                "rateCG", new RealParameter(rates[3].toString()),
                "rateCT", new RealParameter(rates[4].toString()),
                "rateGT", new RealParameter(rates[5].toString()),
                "frequencies", freqs
        );
        */
        nrOfStates = gt16.getStateCount();
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=7)
     * t <- 10
     *
     * rates <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
     * pi <- rep(1/16, 16)
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
     * Q <- sweep(Q, MARGIN=2, pi, `*`)
     *
     * d <- -1 * rowSums(Q)
     * diag(Q) <- d
     * beta <- as.vector(-1 / (pi %*% d))
     * expm(beta * Q * t)
     *
     */
    @Test
    public void testTransitionLong() {
        double t = 10;

        Double[] pi = new Double[nrOfStates];
        Arrays.fill(pi, 1.0 / nrOfStates);

        // rates AC, AG, AT, CG, CT, GT
        Double[] rates = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

        setupModel(pi, rates);

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

        gt16.getTransitionProbabilities(null, t, 0, 1, observed);
        assertArrayEquals(expected, observed, DELTA);
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=7)
     * t <- 0.1
     *
     * rates <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
     * pi <- rep(1/16, 16)
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
     * Q <- sweep(Q, MARGIN=2, pi, `*`)
     *
     * d <- -1 * rowSums(Q)
     * diag(Q) <- d
     * beta <- as.vector(-1 / (pi %*% d))
     * expm(beta * Q * t)
     *
     */
    @Test
    public void testTransitionShort() {
        double t = 0.1;

        Double[] pi = new Double[nrOfStates];
        Arrays.fill(pi, 1.0 / nrOfStates);

        // rates AC, AG, AT, CG, CT, GT
        Double[] rates = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

        setupModel(pi, rates);

        double[] expected = new double[] {
            9.447596e-01, 4.700689e-03, 9.103439e-03, 1.342371e-02, 4.700689e-03, 2.338846e-05, 4.529452e-05, 6.679018e-05, 9.103439e-03, 4.529452e-05, 8.771818e-05, 0.0001293471, 1.342371e-02, 6.679018e-05, 0.0001293471, 0.0001907320,
            4.700689e-03, 9.272311e-01, 1.790894e-02, 2.214672e-02, 2.338846e-05, 4.613475e-03, 8.910663e-05, 1.101919e-04, 4.529452e-05, 8.934539e-03, 1.725655e-04, 0.0002133996, 6.679018e-05, 1.317465e-02, 0.0002544609, 0.0003146738,
            9.103439e-03, 1.790894e-02, 9.185906e-01, 2.638451e-02, 4.529452e-05, 8.910663e-05, 4.570484e-03, 1.312772e-04, 8.771818e-05, 1.725655e-04, 8.851281e-03, 0.0002542337, 1.293471e-04, 2.544609e-04, 0.0130518813, 0.0003748868,
            1.342371e-02, 2.214672e-02, 2.638451e-02, 9.100325e-01, 6.679018e-05, 1.101919e-04, 1.312772e-04, 4.527903e-03, 1.293471e-04, 2.133996e-04, 2.542337e-04, 0.0087688179, 1.907320e-04, 3.146738e-04, 0.0003748868, 0.0129302834,
            4.700689e-03, 2.338846e-05, 4.529452e-05, 6.679018e-05, 9.272311e-01, 4.613475e-03, 8.934539e-03, 1.317465e-02, 1.790894e-02, 8.910663e-05, 1.725655e-04, 0.0002544609, 2.214672e-02, 1.101919e-04, 0.0002133996, 0.0003146738,
            2.338846e-05, 4.613475e-03, 8.910663e-05, 1.101919e-04, 4.613475e-03, 9.100278e-01, 1.757667e-02, 2.173583e-02, 8.910663e-05, 1.757667e-02, 3.394833e-04, 0.0004198151, 1.101919e-04, 2.173583e-02, 0.0004198151, 0.0005191558,
            4.529452e-05, 8.910663e-05, 4.570484e-03, 1.312772e-04, 8.934539e-03, 1.757667e-02, 9.015476e-01, 2.589499e-02, 1.725655e-04, 3.394833e-04, 1.741288e-02, 0.0005001469, 2.133996e-04, 4.198151e-04, 0.0215332791, 0.0006184965,
            6.679018e-05, 1.101919e-04, 1.312772e-04, 4.527903e-03, 1.317465e-02, 2.173583e-02, 2.589499e-02, 8.931483e-01, 2.544609e-04, 4.198151e-04, 5.001469e-04, 0.0172506487, 3.146738e-04, 5.191558e-04, 0.0006184965, 0.0213326642,
            9.103439e-03, 4.529452e-05, 8.771818e-05, 1.293471e-04, 1.790894e-02, 8.910663e-05, 1.725655e-04, 2.544609e-04, 9.185906e-01, 4.570484e-03, 8.851281e-03, 0.0130518813, 2.638451e-02, 1.312772e-04, 0.0002542337, 0.0003748868,
            4.529452e-05, 8.934539e-03, 1.725655e-04, 2.133996e-04, 8.910663e-05, 1.757667e-02, 3.394833e-04, 4.198151e-04, 4.570484e-03, 9.015476e-01, 1.741288e-02, 0.0215332791, 1.312772e-04, 2.589499e-02, 0.0005001469, 0.0006184965,
            8.771818e-05, 1.725655e-04, 8.851281e-03, 2.542337e-04, 1.725655e-04, 3.394833e-04, 1.741288e-02, 5.001469e-04, 8.851281e-03, 1.741288e-02, 8.931464e-01, 0.0256536824, 2.542337e-04, 5.001469e-04, 0.0256536824, 0.0007368461,
            1.293471e-04, 2.133996e-04, 2.542337e-04, 8.768818e-03, 2.544609e-04, 4.198151e-04, 5.001469e-04, 1.725065e-02, 1.305188e-02, 2.153328e-02, 2.565368e-02, 0.8848253781, 3.748868e-04, 6.184965e-04, 0.0007368461, 0.0254146797,
            1.342371e-02, 6.679018e-05, 1.293471e-04, 1.907320e-04, 2.214672e-02, 1.101919e-04, 2.133996e-04, 3.146738e-04, 2.638451e-02, 1.312772e-04, 2.542337e-04, 0.0003748868, 9.100325e-01, 4.527903e-03, 0.0087688179, 0.0129302834,
            6.679018e-05, 1.317465e-02, 2.544609e-04, 3.146738e-04, 1.101919e-04, 2.173583e-02, 4.198151e-04, 5.191558e-04, 1.312772e-04, 2.589499e-02, 5.001469e-04, 0.0006184965, 4.527903e-03, 8.931483e-01, 0.0172506487, 0.0213326642,
            1.293471e-04, 2.544609e-04, 1.305188e-02, 3.748868e-04, 2.133996e-04, 4.198151e-04, 2.153328e-02, 6.184965e-04, 2.542337e-04, 5.001469e-04, 2.565368e-02, 0.0007368461, 8.768818e-03, 1.725065e-02, 0.8848253781, 0.0254146797,
            1.907320e-04, 3.146738e-04, 3.748868e-04, 1.293028e-02, 3.146738e-04, 5.191558e-04, 6.184965e-04, 2.133266e-02, 3.748868e-04, 6.184965e-04, 7.368461e-04, 0.0254146797, 1.293028e-02, 2.133266e-02, 0.0254146797, 0.8765818971
        };

        double[] observed = new double[nrOfStates * nrOfStates];

        gt16.getTransitionProbabilities(null, t, 0, 1, observed);
        assertArrayEquals(expected, observed, DELTA);
    }

}
