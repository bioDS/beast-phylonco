package phylonco.lphy.evolution.readcountmodel;

import lphy.core.logger.LoggerUtils;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.special.Gamma;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ReadCountDataFilter extends DeterministicFunction<ReadCountData> {
    private double fw;
    private double fa;
    private  double ww;
    private double wa;
    private double lambda;
    private double threshold;
    private Value<ReadCountData> readCountData; ;
    public ReadCountDataFilter(
            @ParameterInfo(name = "rc", narrativeName = "read count data", description = "read count data which used to filter variable sites.") Value<ReadCountData> readCountData,
            @ParameterInfo(name = "fw", narrativeName = "expected frequence of alternative nucleotide", description = "the expected frequence of alternative nucleotide in beta-binomial distribution, by default is 0.001", optional = true) Value<Double> fw,
            @ParameterInfo(name = "ww", narrativeName = "shape parameter of homozygous genotype(wild type)", description = "shape parameter of homozygous genotype in beta-binomial distribution, by default is 1.0", optional = true) Value<Double> ww,
            @ParameterInfo(name = "wa", narrativeName = "shape parameter of heterozygous genotype(mutant type)", description = "shape parameter of heterozygous genotype in beta-binomial distribution, by default is 2.0", optional = true) Value<Double> wa,
            @ParameterInfo(name = "lambda", narrativeName = "prior probability of mutation", description = "prior probability of a mutation occurring at a locus, by default is 0.001", optional = true) Value<Double> lambda,
            @ParameterInfo(name = "threshold", narrativeName = "threshold", description = "threshold for determining candidate status, by default is 0.95", optional = true) Value<Double> threshold) {
        this.readCountData = readCountData;
        if (fw != null) {
            this.fw = fw.value();
        } else {
            this.fw = 0.001;
        }

        if (ww != null) {
            this.ww = ww.value();
        } else {
            this.ww = 1.0;
        }

        if (wa != null) {
            this.wa = wa.value();
        } else {
            this.wa = 2.0;
        }

        if (lambda != null) {
            this.lambda = lambda.value();
        } else {
            this.lambda = 0.004;
        }

        if (threshold != null) {
            this.threshold = threshold.value();
        } else {
            this.threshold = 0.05;
        }

        this.fa = 0.5 - this.fw*2/3;
    }


    @GeneratorInfo(name="readCountDataFilter",
            narrativeName = "read count data filter",
            description = "A filter that identify candidate sites from a read counts dataset.")
    @Override
    public Value<ReadCountData> apply() {
        ReadCountData rc = readCountData.value();
        int[] refIndex = rc.getRefIndex();
        int[] orginalIndices = rc.getSitesIndex(); // this is the incides from chrom

        List<Integer> siteIndex = new ArrayList<>();
        for (int i = 0; i < rc.nchar().intValue(); i++) {
            ReadCount[] readCounts = rc.getReadCountsBySite(i);
            int ref = refIndex[i];
            double[] logLikelihoods = logLikelihood(readCounts, ref);
            double[] logProb = logProb(logLikelihoods);
            double[] prob = normalizeLogProbs(logProb);
            if (prob[0] <= threshold) {
                siteIndex.add(i);
            }
        }

        int[] site = new int[siteIndex.size()];
        ReadCount[][] readCountDataMatrix = new ReadCount[rc.getTaxa().ntaxa()][siteIndex.size()];
        for (int i = 0; i < rc.getTaxa().ntaxa(); i++) {
            for (int j = 0; j < siteIndex.size(); j++) {
                readCountDataMatrix[i][j] = rc.getReadCountDataMatrix()[i][siteIndex.get(j)];
                if (j == 0){
                    site[j] = orginalIndices[j];
                }
            }
        }
        ReadCountData filtedRc = new ReadCountData(rc.getTaxa(), readCountDataMatrix, site);

        // print out positions information
        List<Integer> positions = new ArrayList<>();
        for (int i = 0; i < site.length; i++) {
            positions.add(site[i] + 1);
        }
        if (positions.size() == 0) {
            LoggerUtils.log.info("Extract 0 candidate sites");
        } else if (positions.size() == 1){
            LoggerUtils.log.info("Extract " + positions.size() + " candidate sites, the position is " + positions.get(0));
        } else if (positions.size() > 1) {
            LoggerUtils.log.info("Extract " + positions.size() + " candidate sites, the positions are " + positions);
        }

        return new Value<>(null, filtedRc, this);
    }


    public double logProbOfBetaBinomial(int s, int c, double f, double w) {
        double alpha = f * w;
        double beta = w - alpha;
        double logPro = Gamma.logGamma(c+1)+Gamma.logGamma(s+alpha)+Gamma.logGamma(c-s+beta)+Gamma.logGamma(alpha + beta)
                -Gamma.logGamma(s+1)-Gamma.logGamma(c-s+1)-Gamma.logGamma(c+alpha + beta)-Gamma.logGamma(alpha)-Gamma.logGamma(beta);
        return logPro;
    }

    public double logPrior(int m, int k) {
        double logP = 2*logCombinatorialNumber(m, k)-logCombinatorialNumber(2*m, 2*k)-Math.log(2*k-1);
        return logP;
    }

    public double logCombinatorialNumber(int n, int m) {
        double log = Gamma.logGamma(n+1)-Gamma.logGamma(m+1)-Gamma.logGamma(n-m+1);
        return log;
    }

    public double logSumExp(double a, double b) {
        if (Double.isInfinite(a)) return b;
        if (Double.isInfinite(b)) return a;
        double max = Math.max(a, b);
        return max + Math.log(Math.exp(a - max) + Math.exp(b - max));
    }

    public double[] logLikelihood(ReadCount[] readCounts, int ref) {
        int m = readCounts.length;
        double[] logLikelihood = new double[m+1];
        double[] logPwt = new double[readCounts.length];
        double[] logPa = new double[readCounts.length];
        for (int i = 0; i < readCounts.length; i++) {
            int c = readCounts[i].getDepth();
            int r = readCounts[i].getReadCounts()[ref];
            int s = c-r;
            logPwt[i] = logProbOfBetaBinomial(s, c, fw, ww);
            logPa[i] = logProbOfBetaBinomial(s, c, fa, wa);
        }
        for (int k = 0; k <= m; k++) {
            double logP1 = -logCombinatorialNumber(m, k);
            double[][] dp = new double[m + 1][k + 1];
            for (int i = 0; i <= m; i++) {
                Arrays.fill(dp[i], Double.NEGATIVE_INFINITY);
            }
            dp[0][0] = 0.0;
            for (int i = 1; i <= m; i++) {
                for (int j = 0; j <= Math.min(i, k); j++) {
                    double logWild = dp[i - 1][j] + logPwt[i - 1];
                    double logMutant = (j > 0) ? dp[i - 1][j - 1] + logPa[i - 1] : Double.NEGATIVE_INFINITY;
                    dp[i][j] = logSumExp(logWild, logMutant);
                }
            }
            logLikelihood[k] = logP1 + dp[m][k];
        }

        return logLikelihood;
    }

    public double[] logProb(double[] logLikelihood) {
        int m = logLikelihood.length-1;
        double[] logProb = new double[m+1];
        logProb[0] = logLikelihood[0] + Math.log(1-lambda);
        for (int k = 1; k <= m; k++) {
            logProb[k] = logLikelihood[k] + logPrior(m,k) + Math.log(lambda);
        }
        return logProb;
    }

    private double[] normalizeLogProbs(double[] logProbs) {
        // Find the maximum log probability for numerical stability
        double maxLogProb = Double.NEGATIVE_INFINITY;
        for (double logP : logProbs) {
            maxLogProb = Math.max(maxLogProb, logP);
        }
        // Compute exp(logProb - max) for each element and sum them
        double[] expProbs = new double[logProbs.length];
        double sumExp = 0.0;
        for (int i = 0; i < logProbs.length; i++) {
            expProbs[i] = Math.exp(logProbs[i] - maxLogProb);
            sumExp += expProbs[i];
        }
        // Normalize to get actual probabilities
        double[] normalizedProbs = new double[logProbs.length];
        for (int i = 0; i < logProbs.length; i++) {
            normalizedProbs[i] = expProbs[i] / sumExp;
        }
        return normalizedProbs;
    }

}