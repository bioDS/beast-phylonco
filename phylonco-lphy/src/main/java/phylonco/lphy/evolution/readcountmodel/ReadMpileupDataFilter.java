package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.Mpileup;
import lphy.base.evolution.PileupSite;
import lphy.base.function.io.ReadMpileup;
import lphy.base.function.io.ReaderConst;
import lphy.core.io.UserDir;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.special.Gamma;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.*;

import static lphy.base.evolution.PileupSite.translateRead;
import static lphy.base.evolution.datatype.Variant.getCanonicalState;

public class ReadMpileupDataFilter extends ReadMpileup {
    private double fw;
    private double fa;
    private  double ww;
    private double wa;
    private double lambda;
    private double threshold;
    public ReadMpileupDataFilter(@ParameterInfo(name = "fw", narrativeName = "expected sequencing error", description = "the expected sequencing error in beta-binomial distribution, by default is 0.001", optional = true) Value<Double> fw,
                                 @ParameterInfo(name = "ww", narrativeName = "shape parameter of homozygous genotype(wild type)", description = "shape parameter of homozygous genotype in beta-binomial distribution, by default is 1.0", optional = true) Value<Double> ww,
                                 @ParameterInfo(name = "wa", narrativeName = "shape parameter of heterozygous genotype(mutant type)", description = "shape parameter of heterozygous genotype in beta-binomial distribution, by default is 2.0", optional = true) Value<Double> wa,
                                 @ParameterInfo(name = "lambda", narrativeName = "prior probability of mutation", description = "prior probability of a mutation occurring at a locus, by default is 0.0001", optional = true) Value<Double> lambda,
                                 @ParameterInfo(name = "threshold", narrativeName = "threshold", description = "Maximum posterior probability of no mutation (default 0.05, equivalent to requiring mutation posterior > 0.95)", optional = true) Value<Double> threshold,
                                 @ParameterInfo(name = ReaderConst.FILE, description = "the name of mpileup file including path, should have .mpileup or .pileup suffix") Value<String> mpileup,
                                 @ParameterInfo(name = taxaParamName, description = "an array of objects representing taxa names") Value<Object[]> taxaNames) {

        super(mpileup, taxaNames);
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
            this.lambda = 0.0001;
        }

        if (threshold != null) {
            this.threshold = threshold.value();
        } else {
            this.threshold = 0.05;
        }

        this.fa = 0.5 - this.fw*2/3;

    }
    @GeneratorInfo(name = "readMpileupDataFilter", description = "Read in a mpileup file and identify candidate sites. Taxa names are the names used to generate mpileup file.")
    @Override
    public Value<Mpileup> apply() {
        String mpileupFile = getMpileup().value();
        Value names = getTaxaNames();
        String[] taxaNames;
        //Taxa taxaNames = getTaxa().value();

        if (names.value().getClass().isArray()) {
            taxaNames = (String[]) names.value();

        } else throw new IllegalArgumentException(taxaParamName + " must be an array.");


        Mpileup mpileup = readAndFilterMpileup(mpileupFile, taxaNames);

        return new Value<>(null, mpileup, this);
    }

    private Mpileup readAndFilterMpileup(String mpileupFile, String[] taxaNames) {
        Mpileup mpileup;

        Path mpileupPath = Path.of(mpileupFile);

        try (BufferedReader mpileupReader = Files.newBufferedReader(mpileupPath, StandardCharsets.UTF_8)) {

            // Read all cell names into a list

            String line;

            //String[] chroms = new String[mpileupReader.];
            List<String> chroms = new ArrayList<>();
            List<Integer> positions = new ArrayList<>();
            List<Integer> refs = new ArrayList<>();
            List<Map<String, PileupSite.CellPileupData>> pileupData = new ArrayList<>();
            int index = 0;

            while ((line = mpileupReader.readLine()) != null) {
                if (line.isBlank()) continue;

                String[] parts = line.split("\t");
                if (parts.length < 3) {
                    LoggerUtils.log.warning("Skipping malformed line: " + line);
                    continue;
                }

                String chr = parts[0];
                String pos = parts[1];
                String ref = parts[2];


                // Each cell contributes 3 fields: readCount, reads, mapQ
                int expectedDataCols = 3 * taxaNames.length;
                if (parts.length < 3 + expectedDataCols) {
                    LoggerUtils.log.warning("Line has fewer columns than expected for all cells: " + line);
                    continue;
                }

                Map<String, PileupSite.CellPileupData> cellData = new LinkedHashMap<>();
                for (int i = 0; i < taxaNames.length; i++) {
                    int baseIndex = 3 + i * 3;
                    try {
                        int readCount = Integer.parseInt(parts[baseIndex]);
                        String reads = parts[baseIndex + 1];
                        String mapQ = parts[baseIndex + 2];
                        PileupSite.CellPileupData cellPileupData = new PileupSite.CellPileupData(readCount, reads, mapQ);
                        cellData.put(taxaNames[i], cellPileupData);
                    } catch (NumberFormatException e) {
                        LoggerUtils.log.warning("Bad readCount for cell " + taxaNames[i] + " at " + chr + ":" + pos);
                    }
                }
                double[] logLikelihoods = logLikelihood(cellData, taxaNames, getCanonicalState(ref));
                double[] logProb = logProb(logLikelihoods);
                double[] prob = normalizeLogProbs(logProb);
                if (prob[0] <= threshold) {
                    chroms.add(parts[0]);
                    positions.add(Integer.parseInt(parts[1]));
                    refs.add(getCanonicalState(parts[2]));
                    pileupData.add(cellData);
                } else {
                    System.out.println("Not candidate site at Chrom:" + chr + " pos:" + pos);
                }
                index++;
            }
            String[] chrom = chroms.toArray(new String[0]);
            int[] pos = positions.stream().mapToInt(Integer::intValue).toArray();
            int[] ref = refs.stream().mapToInt(Integer::intValue).toArray();
            mpileup = new Mpileup(chrom, pos, ref, pileupData);
            return mpileup;
        } catch (FileNotFoundException | NoSuchFileException e) {
            LoggerUtils.log.severe("File not found: " + e.getMessage() +
                    "\nCurrent working dir = " + UserDir.getUserDir());
        } catch (IOException e) {
            LoggerUtils.logStackTrace(e);
        }

        return null;
    }

    public double[] logLikelihood(Map<String, PileupSite.CellPileupData> cellData, String[] taxaNames, int ref) {
        int m = taxaNames.length;
        double[] logLikelihood = new double[m+1];
        double[] logPwt = new double[m];
        double[] logPa = new double[m];
        for (int i = 0; i < m; i++) {
            PileupSite.CellPileupData cellPileupData = cellData.get(taxaNames[i]);
            String read = translateRead(ref, cellPileupData.reads());
            String[] parts = read.split(":");
            int[] counts = new int[4];
            for (int j = 0; j < counts.length; j++) {
                counts[j] = Integer.parseInt(parts[j].substring(1));
            }
            // Reference canonical nucleotide count.
            int r = counts[ref];
            // Non-reference canonical nucleotide count.
            // Deletions, skips and other non-ACGT observations are excluded.
            int s = 0;
            for (int base = 0; base < counts.length; base++) {
                if (base != ref) {
                    s += counts[base];
                }
            }
            // Effective depth used by the nucleotide model.
            int c = r + s;
            logPwt[i] = logProbOfBetaBinomial(s, c, fw, ww);
            logPa[i] = logProbOfBetaBinomial(s, c, fa, wa);
        }

        double[] prev = new double[m + 1];
        double[] curr = new double[m + 1];
        Arrays.fill(prev, Double.NEGATIVE_INFINITY);
        prev[0] = 0.0;
        for (int i = 1; i <= m; i++) {
            Arrays.fill(curr, Double.NEGATIVE_INFINITY);
            for (int j = 0; j <= i; j++) {
                if (j <= i - 1) {
                    double logWild = prev[j] + logPwt[i - 1];
                    curr[j] = logSumExp(curr[j], logWild);
                }
                if (j > 0) {
                    double logMutant = prev[j - 1] + logPa[i - 1];
                    curr[j] = logSumExp(curr[j], logMutant);
                }
            }
            double[] temp = prev;
            prev = curr;
            curr = temp;
        }
        for (int k = 0; k <= m; k++) {
            double logBinomial = logCombinatorialNumber(m, k);
            logLikelihood[k] = prev[k] - logBinomial;
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
}
