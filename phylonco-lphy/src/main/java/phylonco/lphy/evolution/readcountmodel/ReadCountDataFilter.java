package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.CellPosition;
import lphy.base.evolution.Mpileup;
import lphy.base.evolution.PileupSite;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.special.Gamma;

import java.util.*;

import static lphy.base.evolution.PileupSite.translateRead;
import static phylonco.lphy.evolution.readcountmodel.MpileupToReadCount.translateReads;

public class ReadCountDataFilter extends DeterministicFunction<CellPosition[]> {
    private double fw;
    private double fa;
    private  double ww;
    private double wa;
    private double lambda;
    private double threshold;
    private Value<List<Mpileup>> mpileups;

    public ReadCountDataFilter(
            @ParameterInfo(name = "mpileup", narrativeName = "mpileup input", description = "mpileup input which used to filter candidate sites.") Value<List<Mpileup>> mpileupData,
            @ParameterInfo(name = "fw", narrativeName = "expected frequence of alternative nucleotide", description = "the expected frequence of alternative nucleotide in beta-binomial distribution, by default is 0.001", optional = true) Value<Double> fw,
            @ParameterInfo(name = "ww", narrativeName = "shape parameter of homozygous genotype(wild type)", description = "shape parameter of homozygous genotype in beta-binomial distribution, by default is 1.0", optional = true) Value<Double> ww,
            @ParameterInfo(name = "wa", narrativeName = "shape parameter of heterozygous genotype(mutant type)", description = "shape parameter of heterozygous genotype in beta-binomial distribution, by default is 2.0", optional = true) Value<Double> wa,
            @ParameterInfo(name = "lambda", narrativeName = "prior probability of mutation", description = "prior probability of a mutation occurring at a locus, by default is 0.001", optional = true) Value<Double> lambda,
            @ParameterInfo(name = "threshold", narrativeName = "threshold", description = "threshold for determining candidate status, by default is 0.95", optional = true) Value<Double> threshold) {
        if (mpileupData == null) {
            throw new NullPointerException("mpileup data is null");
        }
        this.mpileups = mpileupData;
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
            narrativeName = "read count data filter", examples = {"mpileupToReadCount.lphy"},
            description = "A filter that identify candidate sites from a read counts dataset.")
    @Override
    public Value<CellPosition[]> apply() {
        List<Mpileup> mpileupData = getMpileups().value();

        Map<Integer, List<ReadCount>> positionToReadCounts = new LinkedHashMap<>();
        Map<Integer, List<CellPosition>> positionToCellPositions = new LinkedHashMap<>();
        Map<Integer, Integer> positionToRef = new LinkedHashMap<>();

        for (Mpileup mp : mpileupData) {
            int pos = mp.getPosition();
            int ref = mp.getRef();
            String chromName = mp.getChromName();

            positionToRef.put(pos, ref);

            List<ReadCount> rcList = positionToReadCounts.computeIfAbsent(pos, k -> new ArrayList<>());
            List<CellPosition> cpList = positionToCellPositions.computeIfAbsent(pos, k -> new ArrayList<>());

            for (Map.Entry<String, PileupSite.CellPileupData> entry : mp.getPileupData().entrySet()) {
                String cellName = entry.getKey();
                PileupSite.CellPileupData pileup = entry.getValue();

                ReadCount rc = translateReads(ref, pileup);
                rcList.add(rc);

                CellPosition cp = new CellPosition(chromName, cellName, pos);
                cpList.add(cp);
            }
        }

        List<Integer> siteIndex = new ArrayList<>();
        for (int pos : positionToReadCounts.keySet()) {
            List<ReadCount> reads = positionToReadCounts.get(pos);
            int ref = positionToRef.get(pos);
            double[] logLikelihoods = logLikelihood(reads, ref);
            double[] logProb = logProb(logLikelihoods);
            double[] prob = normalizeLogProbs(logProb);
            if (prob[0] <= threshold) {
                siteIndex.add(pos);
            }
        }

        List<CellPosition> cellPositions = new ArrayList<>();
        for (Integer site : siteIndex) {
            cellPositions.addAll(positionToCellPositions.get(site));
        }

        if (siteIndex.size() == 0) {
            LoggerUtils.log.info("Extract 0 candidate sites");
        } else if (siteIndex.size() == 1){
            LoggerUtils.log.info("Extract " + siteIndex.size() + " candidate sites, the position is " + siteIndex.get(0));
        } else if (siteIndex.size() > 1) {
            LoggerUtils.log.info("Extract " + siteIndex.size() + " candidate sites, the positions are " + siteIndex);
        }

        return new Value<>("", cellPositions.toArray(new CellPosition[0]), this);
    }

    public Value<List<Mpileup>> getMpileups() {
        return getParams().get("mpileup");
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

    public double[] logLikelihood(List<ReadCount> readCounts, int ref) {
        int m = readCounts.size();
        double[] logLikelihood = new double[m+1];
        double[] logPwt = new double[readCounts.size()];
        double[] logPa = new double[readCounts.size()];
        for (int i = 0; i < readCounts.size(); i++) {
            int c = readCounts.get(i).getDepth();
            int r = readCounts.get(i).getReadCounts()[ref];
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
