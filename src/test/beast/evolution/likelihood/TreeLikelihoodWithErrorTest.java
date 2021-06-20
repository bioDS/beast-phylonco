package test.beast.evolution.likelihood;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.Binary;
import beast.evolution.datatype.Nucleotide;
import beast.evolution.datatype.NucleotideDiploid16;
import beast.evolution.errormodel.BinaryErrorModel;
import beast.evolution.errormodel.ErrorModelBase;
import beast.evolution.errormodel.GT16ErrorModel;
import beast.evolution.likelihood.TreeLikelihoodWithError;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GT16;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.substitutionmodel.BinarySubstitutionModel;
import beast.util.TreeParser;
import org.junit.Test;

import java.util.Arrays;

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
     *   c(1.0 - ep, ep/3, ep/3, ep/3,
     *   ep/3, 1.0 - ep, ep/3, ep/3,
     *   ep/3, ep/3, 1.0 - ep, ep/3,
     *   ep/3, ep/3, ep/3, 1.0 - ep), nrow=4, byrow=T)
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
                "dataType", "nucleotide"
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

        Nucleotide datatype = new Nucleotide();

        ErrorModelBase errorModel = new ErrorModelBase();
        errorModel.initByName("epsilon", "0.1", "datatype", datatype);
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
        double expectedLogP = -2.3063595712034233;
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
                "dataType", "binary"
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

        Binary datatype = new Binary();

        BinaryErrorModel errorModel = new BinaryErrorModel();
        errorModel.initByName("alpha", alpha, "beta", beta, "datatype", datatype);
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

    private double calculateLikelihoodGT16(String seq, String epsilon, String delta) {
        Alignment data = new Alignment();
        Sequence seqA = new Sequence("a", seq.substring(0, 1));
        Sequence seqB = new Sequence("b", seq.substring(1));
        data.initByName(
                "sequence", seqA,
                "sequence", seqB,
                "dataType", "nucleotideDiploid16"
        );

        TreeParser tree = new TreeParser();
        tree.initByName(
                "taxa", data,
                "newick", "(a: 0.5, b: 0.5);",
                "IsLabelledNewick", true
        );

        Double[] pi = new Double[16];
        Arrays.fill(pi, 1.0 / 16);
        RealParameter f = new RealParameter(pi);
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f, "estimate", false);
        freqs.initAndValidate();

        Double[] rates = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        RealParameter nucRates = new RealParameter(rates);
        nucRates.setInputValue("keys", "AC AG AT CG CT GT");
        nucRates.initAndValidate();

        GT16 subsModel = new GT16();
        subsModel.initByName(
                "nucRates", nucRates,
                "frequencies", freqs
        );
        subsModel.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", subsModel);
        siteModel.initAndValidate();

        NucleotideDiploid16 datatype = new NucleotideDiploid16();

        GT16ErrorModel errorModel = new GT16ErrorModel();
        errorModel.initByName("epsilon", epsilon, "delta", delta, "datatype", datatype);
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

    /***
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=7)
     * t <- 0.5
     *
     * # substitution model
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
     * P <- expm(beta * Q * t)
     *
     * # error model
     * ep <- 0.1
     * dt <- 0.2
     *
     * a <- 1 - ep + (1/2.0) * dt * ep
     * b <- (1 - dt) * (1/6.0) * ep
     * c <- (1/6.0) * dt * ep
     * d <- (1/2.0) * dt + (1/6.0) * ep - (1/3.0) * dt * ep
     * e <- (1/6.0) * dt * ep
     * f <- (1 - dt) * (1/6.0) * ep
     * g <- (1 - dt) * (1 - ep)
     *
     * err <- matrix(c(
     * #	AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
     * 	a, b, b, b, b, c, 0, 0, b, 0, c, 0, b, 0, 0, c, # AA
     * 	d, g, f, f, 0, d, 0, 0, 0, f, e, 0, 0, f, 0, e, # AC
     * 	d, f, g, f, 0, e, f, 0, 0, 0, d, 0, 0, 0, f, e, # AG
     * 	d, f, f, g, 0, e, 0, f, 0, 0, e, f, 0, 0, 0, d, # AT
     * 	d, 0, 0, 0, g, d, f, f, f, 0, e, 0, f, 0, 0, e, # CA
     * 	c, b, 0, 0, b, a, b, b, 0, b, c, 0, 0, b, 0, c, # CC
     * 	e, 0, f, 0, f, d, g, f, 0, 0, d, 0, 0, 0, f, e, # CG
     * 	e, 0, 0, f, f, d, f, g, 0, 0, e, f, 0, 0, 0, d, # CT
     * 	d, 0, 0, 0, f, e, 0, 0, g, f, d, f, f, 0, 0, e, # GA
     * 	e, f, 0, 0, 0, d, 0, 0, f, g, d, f, 0, f, 0, e, # GC
     * 	c, 0, b, 0, 0, c, b, 0, b, b, a, b, 0, 0, b, c, # GG
     * 	e, 0, 0, f, 0, e, 0, f, f, f, d, g, 0, 0, 0, d, # GT
     * 	d, 0, 0, 0, f, e, 0, 0, f, 0, e, 0, g, f, f, d, # TA
     * 	e, f, 0, 0, 0, d, 0, 0, 0, f, e, 0, f, g, f, d, # TC
     * 	e, 0, f, 0, 0, e, f, 0, 0, 0, d, 0, f, f, g, d, # TG
     * 	c, 0, 0, b, 0, c, 0, b, 0, 0, c, b, b, b, b, a  # TT
     * ), nrow=16, byrow=F)
     * p1 <- P %*% err[1,]
     * prob <- 0.0
     * for (i in 1:16) {
     * 	prob <- pi[i] * (p1[i] ** 2) + prob
     * }
     * log(prob)
     */
    @Test
    public void testGT16ErrorLikelihoodCase0() {
        double logP = calculateLikelihoodGT16("00", "0.1", "0.2");
        double expectedLogP = -3.2683402019565975;
        assertEquals(expectedLogP, logP, DELTA);
    }

    /***
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=7)
     * t <- 0.5
     *
     * # substitution model
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
     * P <- expm(beta * Q * t)
     *
     * # error model
     * ep <- 0.1
     * dt <- 0.2
     *
     * a <- 1 - ep + (1/2.0) * dt * ep
     * b <- (1 - dt) * (1/6.0) * ep
     * c <- (1/6.0) * dt * ep
     * d <- (1/2.0) * dt + (1/6.0) * ep - (1/3.0) * dt * ep
     * e <- (1/6.0) * dt * ep
     * f <- (1 - dt) * (1/6.0) * ep
     * g <- (1 - dt) * (1 - ep)
     *
     * err <- matrix(c(
     * #	AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
     * 	a, b, b, b, b, c, 0, 0, b, 0, c, 0, b, 0, 0, c, # AA
     * 	d, g, f, f, 0, d, 0, 0, 0, f, e, 0, 0, f, 0, e, # AC
     * 	d, f, g, f, 0, e, f, 0, 0, 0, d, 0, 0, 0, f, e, # AG
     * 	d, f, f, g, 0, e, 0, f, 0, 0, e, f, 0, 0, 0, d, # AT
     * 	d, 0, 0, 0, g, d, f, f, f, 0, e, 0, f, 0, 0, e, # CA
     * 	c, b, 0, 0, b, a, b, b, 0, b, c, 0, 0, b, 0, c, # CC
     * 	e, 0, f, 0, f, d, g, f, 0, 0, d, 0, 0, 0, f, e, # CG
     * 	e, 0, 0, f, f, d, f, g, 0, 0, e, f, 0, 0, 0, d, # CT
     * 	d, 0, 0, 0, f, e, 0, 0, g, f, d, f, f, 0, 0, e, # GA
     * 	e, f, 0, 0, 0, d, 0, 0, f, g, d, f, 0, f, 0, e, # GC
     * 	c, 0, b, 0, 0, c, b, 0, b, b, a, b, 0, 0, b, c, # GG
     * 	e, 0, 0, f, 0, e, 0, f, f, f, d, g, 0, 0, 0, d, # GT
     * 	d, 0, 0, 0, f, e, 0, 0, f, 0, e, 0, g, f, f, d, # TA
     * 	e, f, 0, 0, 0, d, 0, 0, 0, f, e, 0, f, g, f, d, # TC
     * 	e, 0, f, 0, 0, e, f, 0, 0, 0, d, 0, f, f, g, d, # TG
     * 	c, 0, 0, b, 0, c, 0, b, 0, 0, c, b, b, b, b, a  # TT
     * ), nrow=16, byrow=F)
     * p1 <- P %*% err[1,]
     * p2 <- P %*% err[2,]
     * prob <- 0.0
     * for (i in 1:16) {
     * 	prob <- pi[i] * (p1[i] * p2[i]) + prob
     * }
     * log(prob)
     */
    @Test
    public void testGT16ErrorLikelihoodCase1() {
        double logP = calculateLikelihoodGT16("01", "0.1", "0.2");
        double expectedLogP = -5.1071258693509041;
        assertEquals(expectedLogP, logP, DELTA);
    }
}
