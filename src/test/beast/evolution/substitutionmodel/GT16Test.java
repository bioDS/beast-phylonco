package test.beast.evolution.substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GT16;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

public class GT16Test {

    private static double DELTA = 1e-10;

    private GT16 model;
    private int nrOfStates;

    public void setupModel(Double[] pi, Double[] rates) {
        model = new GT16();

        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", new RealParameter(pi), "estimate", false);

        model.initByName(
                "frequencies", frequencies,
                "rateAC", new RealParameter(rates[0].toString()),
                "rateAG", new RealParameter(rates[1].toString()),
                "rateAT", new RealParameter(rates[2].toString()),
                "rateCG", new RealParameter(rates[3].toString()),
                "rateCT", new RealParameter(rates[4].toString()),
                "rateGT", new RealParameter(rates[5].toString())
        );
        nrOfStates = model.getStateCount();
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=12)
     * t <- 10
     *
     * rates <- c(1.0, 2.0, 1.0, 1.0, 2.0, 1.0)
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
     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     * ), nrow=2, byrow=T)
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

        Double equalProb = 1.0 / nrOfStates;

        Double[] pi = new Double[nrOfStates];
        Arrays.fill(pi, equalProb);

        // rates AC, AG, AT, CG, CT, GT
        Double[] rates = {1.0, 2.0, 1.0, 1.0, 2.0, 1.0};

        setupModel(pi, rates);

        double[] expected = new double[] {

        };

        double[] observed = new double[nrOfStates * nrOfStates];

        model.getTransitionProbabilities(null, t, 0, 1, observed);
        assertArrayEquals(expected, observed, DELTA);
    }

}
