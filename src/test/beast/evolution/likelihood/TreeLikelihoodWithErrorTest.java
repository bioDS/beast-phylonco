package test.beast.evolution.likelihood;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.BinaryWithError;
import beast.evolution.datatype.NucleotideWithError;
import beast.evolution.likelihood.TreeLikelihoodWithError;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.substitutionmodel.BinarySubstitutionModel;
import beast.util.TreeParser;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class TreeLikelihoodWithErrorTest {

    private static double DELTA = 1e-10;

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=20)
     * mu <- 1.0
     * t <- 0.5
     * delta <- (4.0 / 3.0) * t
     * p <- 0.25 * (1.0 + 3.0 * exp(-delta * mu))
     * q <- 0.25 * (1.0 - exp(-delta * mu))
     * P <- matrix(
     *   c(p, q, q, q,
     *   q, p, q, q,
     *   q, q, p, q,
     *   q, q, q, p), nrow=4, byrow=T)
     *
     * ep <- 0.1
     * err <- matrix(
     *   c(1.0 - 3.0 * ep, ep, ep, ep,
     *   ep, 1.0 - 3.0 * ep, ep, ep,
     *   ep, ep, 1.0 - 3.0 * ep, ep,
     *   ep, ep, ep, 1.0 - 3.0 * ep), nrow=4, byrow=T)
     * p1 <- P %*% err[1,]
     * prob <- 0.25 * (p1[1] ** 2 + p1[2] ** 2 + p1[3] ** 2 + p1[4] ** 2)
     * log(prob)
     */
    @Test
    public void testJCLikelihoodSmallWithError() {
        Alignment data = new Alignment();
        Sequence seqA = new Sequence("a", "A");
        Sequence seqB = new Sequence("b", "A");
        data.initByName(
                "sequence", seqA,
                "sequence", seqB,
                "dataType", "nucleotideWithError"
        );

        TreeParser tree = new TreeParser();
        tree.initByName(
                "taxa", data,
                "newick", "(a: 0.5, b: 0.5);",
                "IsLabelledNewick", true
        );

        JukesCantor subsModel = new JukesCantor();
        subsModel.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", subsModel);
        siteModel.initAndValidate();

        NucleotideWithError errorModel = new NucleotideWithError();
        errorModel.initByName("epsilon", "0.1");
        errorModel.initAndValidate();

        TreeLikelihoodWithError likelihood = new TreeLikelihoodWithError();
        likelihood.initByName(
                "data", data,
                "tree", tree,
                "siteModel", siteModel,
                "useAmbiguities", true,
                "useTipLikelihoods", true,
                "errorModel", errorModel);

        double logP = likelihood.calculateLogP();
        double expectedLogP = -2.5220752408362181;
        assertEquals(expectedLogP, logP, DELTA);

        System.out.println("seq A: " + data.getSequenceAsString("a"));
        System.out.println("seq B: " + data.getSequenceAsString("b"));
        System.out.println("likelihood: " + logP);
    }

    private double calculateLikelihoodBinary(String seq, String alpha, String beta) {
        Alignment data = new Alignment();
        Sequence seqA = new Sequence("a", seq.substring(0, 1));
        Sequence seqB = new Sequence("b", seq.substring(1));
        data.initByName(
                "sequence", seqA,
                "sequence", seqB,
                "dataType", "binaryWithError"
        );

        TreeParser tree = new TreeParser();
        tree.initByName(
                "taxa", data,
                "newick", "(a: 0.5, b: 0.5);",
                "IsLabelledNewick", true
        );

        BinarySubstitutionModel subsModel = new BinarySubstitutionModel();
        subsModel.initByName("lambda", "2.0");
        subsModel.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", subsModel);
        siteModel.initAndValidate();

        BinaryWithError errorModel = new BinaryWithError();
        errorModel.initByName("alpha", alpha, "beta", beta);
        errorModel.initAndValidate();

        TreeLikelihoodWithError likelihood = new TreeLikelihoodWithError();
        likelihood.initByName(
                "data", data,
                "tree", tree,
                "siteModel", siteModel,
                "useAmbiguities", true,
                "useTipLikelihoods", true,
                "errorModel", errorModel);

        return likelihood.calculateLogP();
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=20)
     * lambda <- 2.0
     * t <- 0.5
     * Q <- matrix(c(
     *     -1, 1,
     *     lambda, -lambda
     * ), nrow=2, byrow=T)
     * pi0 <- lambda / (lambda + 1)
     * pi1 <- 1 / (lambda + 1)
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * P <- expm(beta * Q * t)
     *
     * alpha <- 0
     * beta <- 0
     * err <- matrix(
     *   c(1 - alpha, beta,
     *   alpha, 1 - beta), nrow=2, byrow=T)
     * p1 <- P %*% err[1,]
     * prob <- pi0 * (p1[1] ** 2) + pi1 * (p1[2] ** 2)
     * log(prob)
     */
    @Test
    public void testBinaryLikelihoodSmallNoError() {
        double logP = calculateLikelihoodBinary("00", "0.0", "0.0");
        double expectedLogP = -0.7595722922504291;
        assertEquals(expectedLogP, logP, DELTA);
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=20)
     * lambda <- 2.0
     * t <- 0.5
     * Q <- matrix(c(
     *     -1, 1,
     *     lambda, -lambda
     * ), nrow=2, byrow=T)
     * pi0 <- lambda / (lambda + 1)
     * pi1 <- 1 / (lambda + 1)
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * P <- expm(beta * Q * t)
     *
     * alpha <- 0.1
     * beta <- 0.2
     * err <- matrix(
     *   c(1 - alpha, beta,
     *   alpha, 1 - beta), nrow=2, byrow=T)
     * p1 <- P %*% err[1,]
     * prob <- pi0 * (p1[1] ** 2) + pi1 * (p1[2] ** 2)
     * log(prob)
     */
    @Test
    public void testBinaryLikelihoodSmallWithErrorCase0() {
        double logP = calculateLikelihoodBinary("00", "0.1", "0.2");
        double expectedLogP = -0.78543518416993563;
        assertEquals(expectedLogP, logP, DELTA);
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=20)
     * lambda <- 2.0
     * t <- 0.5
     * Q <- matrix(c(
     *     -1, 1,
     *     lambda, -lambda
     * ), nrow=2, byrow=T)
     * pi0 <- lambda / (lambda + 1)
     * pi1 <- 1 / (lambda + 1)
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * P <- expm(beta * Q * t)
     *
     * alpha <- 0.1
     * beta <- 0.2
     * err <- matrix(
     *   c(1 - alpha, beta,
     *   alpha, 1 - beta), nrow=2, byrow=T)
     * p1 <- P %*% err[2,]
     * prob <- pi0 * (p1[1] ** 2) + pi1 * (p1[2] ** 2)
     * log(prob)
     */
    @Test
    public void testBinaryLikelihoodSmallWithErrorCase1() {
        double logP = calculateLikelihoodBinary("11", "0.1", "0.2");
        double expectedLogP = -2.0989268283365146;
        assertEquals(expectedLogP, logP, DELTA);
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=20)
     * lambda <- 2.0
     * t <- 0.5
     * Q <- matrix(c(
     *     -1, 1,
     *     lambda, -lambda
     * ), nrow=2, byrow=T)
     * pi0 <- lambda / (lambda + 1)
     * pi1 <- 1 / (lambda + 1)
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * P <- expm(beta * Q * t)
     *
     * alpha <- 0.1
     * beta <- 0.2
     * err <- matrix(
     *   c(1 - alpha, beta,
     *   alpha, 1 - beta), nrow=2, byrow=T)
     * p1 <- P %*% err[1,]
     * p2 <- P %*% err[2,]
     * prob <- pi0 * (p1[1] * p2[1]) + pi1 * (p1[2] * p2[2])
     * log(prob)
     */
    @Test
    public void testBinaryLikelihoodSmallWithErrorCase2() {
        double logP = calculateLikelihoodBinary("01", "0.1", "0.2");
        double expectedLogP = -1.5571044248279775;
        assertEquals(expectedLogP, logP, DELTA);
    }

    /**
    /* tests total probability of all data permutations sums to 1.0
    */
    @Test
    public void testBinaryLikelihoodSmallTotalProbability() {
        double logP1 = calculateLikelihoodBinary("00", "0.1", "0.2");
        double logP2 = calculateLikelihoodBinary("01", "0.1", "0.2");
        double logP3 = calculateLikelihoodBinary("10", "0.1", "0.2");
        double logP4 = calculateLikelihoodBinary("11", "0.1", "0.2");
        double probSum = Math.exp(logP1) + Math.exp(logP2) + Math.exp(logP3) + Math.exp(logP4);
        assertEquals(1.0, probSum, DELTA);
    }
}
