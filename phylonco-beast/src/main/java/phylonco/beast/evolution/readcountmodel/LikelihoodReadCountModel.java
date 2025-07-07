package phylonco.beast.evolution.readcountmodel;



import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import mutablealignment.MutableAlignment;
import org.apache.commons.math3.special.Gamma;
import phylonco.beast.evolution.datatype.ReadCount;

import java.util.HashMap;
import java.util.List;
import java.util.Random;


public class LikelihoodReadCountModel extends Distribution {

    public Input<Alignment> alignmentInput = new Input<>("alignment", "alignment");
    public Input<ReadCount> readCountInput = new Input<>("readCount", "nucleotide read counts");

    // epsilon, allelic dropout, ... parameters
    public Input<RealParameter> epsilonInput = new Input<>("epsilon", "sequencing error");
    public Input<RealParameter> deltaInput = new Input<>("delta", "allelic dropout probability");
    public Input<RealParameter> tInput = new Input<>("t", "mean of allelic coverage");
    public Input<RealParameter> vInput = new Input<>("v", "variance of allelic coverage");
    public Input<RealParameter> sInput = new Input<>("s", "size factor of cell");
    public Input<RealParameter> w1Input = new Input<>("w1", "homozygous genotype overdispersion parameter of Dirichlet multinomial distribution");
    public Input<RealParameter> w2Input = new Input<>("w2", "heterogeneous genotype overdispersion parameter of Dirichlet multinomial distribution");

    // other parameters of read count model

    private RealParameter epsilon;
    private RealParameter delta;
    private RealParameter t;
    private RealParameter v;
    private RealParameter s;
    private RealParameter w1;
    private RealParameter w2;
    private Alignment alignment;
    private ReadCount readCount;
    private double[] negp1, negp2, negr1, negr2;
    private double[] w;
    private Double[][] propensities;
    private int[][] coverages;
    private final double[] alpha = new double[]{1.0, 2.0};
    private double [] currentLogPi, storedLogPi;
    private HashMap<Integer, Double> logGammaCache;
    private HashMap<Integer, Double> logCache;
    private double pLog;
    private double p1Log;
    private double[] wLogGamma = new double[2];
    private double[] deltaLog = new double[2];
    private final double log0_5 = Math.log(0.5);



    @Override
    public List<String> getArguments() {
        return List.of();
    }

    @Override
    public List<String> getConditions() {
        return List.of();
    }

    @Override
    public void sample(State state, Random random) {

    }

    @Override
    public void initAndValidate() {
        // checking parameters correct
        // get parameters
        epsilon = epsilonInput.get();
        delta = deltaInput.get();
        t = tInput.get();
        v = vInput.get();
        s = sInput.get();
        w1 = w1Input.get();
        w2 = w2Input.get();
        alignment = alignmentInput.get();
        readCount = readCountInput.get();
        negp1 = new double[s.getDimension()];
        negp2 = new double[s.getDimension()];
        negr1 = new double[s.getDimension()];
        negr2 = new double[s.getDimension()];
        coverages = new int[alignment.getTaxonCount()][alignment.getSiteCount()];
        for (int i = 0; i < alignment.getTaxonCount(); i++) {
            for (int j = 0; j < alignment.getSiteCount(); j++) {
                for (int k = 0; k < 4; k++) {
                    coverages[i][j] += readCount.getReadCounts(i,j)[k];
                }
            }
        }


        currentLogPi = new double[alignment.getTaxonCount()];
        storedLogPi =  new double[alignment.getTaxonCount()];
        logGammaCache = new HashMap<>();
        logCache = new HashMap<>();
    }

    // calculate propensities matrix of dirichlet multinomial distribution(read count model)
    // and params of negative binomial distribution(coverage model)
    private void initialize() {
        double mean1;
        double mean2;
        double variance1;
        double variance2;
        Double eps = epsilon.getValue();
        Double tv = t.getValue();
        Double vv = v.getValue();
        Double[] sv = s.getValues();

        for (int i = 0; i < s.getDimension(); i++) {
            mean1 = alpha[0] * tv * sv[i];
            mean2 = alpha[1] * tv * sv[i];
            variance1 = mean1 + Math.pow(alpha[0], 2) * vv * Math.pow(sv[i], 2);
            variance2 = mean2 + Math.pow(alpha[1], 2) * vv * Math.pow(sv[i], 2);
            negp1[i] = mean1 / variance1;
            negp2[i] = mean2 / variance2;
            negr1[i] = Math.pow(mean1, 2) / (variance1 - mean1);
            negr2[i] = Math.pow(mean2, 2) / (variance2 - mean2);
        }
        w = new double[]{w1.getValue(), w2.getValue()};
        propensities = new Double[][]{
                {(1 - eps), (eps/3), (eps/3), (eps/3)},             // AA or A_ 0
                {(0.5 - eps/6), (0.5 - eps/6), (eps/6), (eps/6)},   // AC or CA 1
                {(0.5 - eps/6), (eps/6), (0.5 - eps/6), (eps/6)},   // AG or GA 2
                {(0.5 - eps/6), (eps/6), (eps/6), (0.5 - eps/6)},   // AT or TA 3
                {(eps/3), (1 - eps), (eps/3), (eps/3)},             // CC or C_ 4
                {(eps/6), (0.5 - eps/6), (0.5 - eps/6), (eps/6)},   // CG or GC 5
                {(eps/6), (0.5 - eps/6), (eps/6), (0.5 - eps/6)},   // CT or TC 6
                {(eps/3), (eps/3), (1 - eps), (eps/3)},             // GG or G_ 7
                {(eps/6), (eps/6), (0.5 - eps/6), (0.5 - eps/6)},   // GT or TG 8
                {(eps/3), (eps/3), (eps/3), (1 - eps)},             // TT or T_ 9
        };

        wLogGamma[0] = Gamma.logGamma(w[0]);
        wLogGamma[1] = Gamma.logGamma(w[1]);
        deltaLog[0] = Math.log(delta.getValue());
        deltaLog[1] = Math.log(1-delta.getValue());

    }

    public void calculateLogPLeaf(Node node, int[] states) {

    }

    //Calculate the log likelihood of read count model by summarizing the log likelihood at each site
    @Override
    public double calculateLogP() {
        initialize();
        if (alignment instanceof MutableAlignment a) {
            return calculateLogP(a);
        }
        //logP = 0;
        for (int i = 0; i < alignment.getTaxonCount(); i++) {
            double logPi = 0;
            for (int j = 0; j < alignment.getSiteCount(); j++) {///？
                // dirichlet multinomial pmf
                int patternIndex = alignment.getPatternIndex(j);
                int genotypeState = alignment.getPattern(i, patternIndex);
                int[] readCountNumbers = readCount.getReadCounts(i, j);
                logPi += logLiklihoodRC(genotypeState, readCountNumbers, coverages[i][j], i);
            }
            currentLogPi[i] = logPi;
        }
        logP = 0;
        for (double d : currentLogPi) {
            logP += d;
        }
        return logP;
    }

    private double calculateLogP(MutableAlignment mutableAlignment) {
        /** update currentLogPi only for sequences that changed **/
        if (mutableAlignment.getDirtySequenceIndices().length != 0) {
            for (int i : mutableAlignment.getDirtySequenceIndices()) {
                double logPi = 0;
                for (int j = 0; j < mutableAlignment.getSiteCount(); j++) {///？
                    // dirichlet multinomial pmf
                    int patternIndex = mutableAlignment.getPatternIndex(j);
                    int genotypeState = mutableAlignment.getPattern(i, patternIndex);
                    int[] readCountNumbers = readCountInput.get().getReadCounts(i, j);
                    logPi += logLiklihoodRC(genotypeState, readCountNumbers, coverages[i][j], i);
                }
                currentLogPi[i] = logPi;
            }

            /** sum over all sequence contributions **/
            logP = 0;
            for (double d : currentLogPi) {
                logP += d;
            }
            return logP;
        } else {
            for (int i = 0; i < alignment.getTaxonCount(); i++) {
                double logPi = 0;
                for (int j = 0; j < alignment.getSiteCount(); j++) {///？
                    // dirichlet multinomial pmf
                    int patternIndex = alignment.getPatternIndex(j);
                    int genotypeState = alignment.getPattern(i, patternIndex);
                    int[] readCountNumbers = readCount.getReadCounts(i, j);
                    logPi += logLiklihoodRC(genotypeState, readCountNumbers, coverages[i][j], i);
                }
                currentLogPi[i] = logPi;
            }
            logP = 0;
            for (double d : currentLogPi) {
                logP += d;
            }
            return logP;
        }

    }

    public double[] sequenceLogLikelihood(int taxonIndex, int[] genotypeSequence) {
        double[] taxonLogP = new double[genotypeSequence.length];
        if (genotypeSequence.length != alignment.getSiteCount()) {
            throw new RuntimeException("genotypeSequence.length != alignment.getSiteCount()");
        }
        for (int j = 0; j < genotypeSequence.length; j++) {
            // dirichlet multinomial pmf
            int genotypeState = genotypeSequence[j];
            int[] readCountNumbers = readCount.getReadCounts(taxonIndex, j);
            taxonLogP[j] = logLiklihoodRC(genotypeState, readCountNumbers, coverages[taxonIndex][j], taxonIndex);
        }
        return taxonLogP;
    }


    // calculate probability of read counts given genotype
    // genotypeState represents genotype alignment
    public double logLiklihoodRC(int genotypeState, int[] readCountNumbers, int coverage, int taxonIndex) {

        int[] indices = getGenotypeIndices(genotypeState);

        double logLikelihood;
        double logLikelihoodDirichletMDDiploid;
        double logLikelihoodDirichletMDHaploid0;
        double logLikelihoodDirichletMDHaploid1;
        double logCoverageLikelihoodDiploid;
        double logCoverageLikelihoodHaploid;
        double part0;
        double part1;
        double part2;
        double max;


        if (homozygous(genotypeState)) {
            logLikelihoodDirichletMDDiploid = logLikelihoodDirichletMD(w[0], coverage, propensities[indices[0]], readCountNumbers);
            logCoverageLikelihoodDiploid = logCoverageLikelihood(coverage,negp2[taxonIndex], negr2[taxonIndex]);
            logLikelihoodDirichletMDHaploid0 = logLikelihoodDirichletMD(w[0], coverage, propensities[indices[1]], readCountNumbers);
            logCoverageLikelihoodHaploid = logCoverageLikelihood(coverage,negp1[taxonIndex], negr1[taxonIndex]);
            part0 = logLikelihoodDirichletMDDiploid + logCoverageLikelihoodDiploid + deltaLog[1];
            part1 = logLikelihoodDirichletMDHaploid0 + logCoverageLikelihoodHaploid + deltaLog[0];
            max = Math.max(part0, part1);
            if (part0 == max){
                logLikelihood = part0 + Math.log(1 + Math.exp(part1 - part0));
            }else {
                logLikelihood = part1 + Math.log(1 + Math.exp(part0 - part1));
            }
        } else {
            logLikelihoodDirichletMDDiploid = logLikelihoodDirichletMD(w[1], coverage, propensities[indices[0]], readCountNumbers);
            logCoverageLikelihoodDiploid = logCoverageLikelihood(coverage, negp2[taxonIndex], negr2[taxonIndex]);
            logLikelihoodDirichletMDHaploid0 = logLikelihoodDirichletMD(w[0], coverage, propensities[indices[1]], readCountNumbers);
            logCoverageLikelihoodHaploid = logCoverageLikelihood(coverage, negp1[taxonIndex], negr1[taxonIndex]);
            logLikelihoodDirichletMDHaploid1 = logLikelihoodDirichletMD(w[0], coverage, propensities[indices[2]], readCountNumbers);
            part0 = logLikelihoodDirichletMDDiploid + logCoverageLikelihoodDiploid + deltaLog[1];
            part1 = log0_5 + logLikelihoodDirichletMDHaploid0 + logCoverageLikelihoodHaploid + deltaLog[0];
            part2 = log0_5 + logLikelihoodDirichletMDHaploid1 + logCoverageLikelihoodHaploid + deltaLog[0];
            max = Math.max(part0, Math.max(part1, part2));
            if (part0 == max){
                logLikelihood = part0 + Math.log(1 + Math.exp(part1 - part0)) + Math.log(1 + Math.exp(part2 - part0));
            }else if (part1 == max){
                logLikelihood = part1 + Math.log(1 + Math.exp(part0 - part1)) + Math.log(1 + Math.exp(part2 - part1));
            }else {
                logLikelihood = part2 + Math.log(1 + Math.exp(part0 - part2)) + Math.log(1 + Math.exp(part1 - part2));
            }
        }
        return logLikelihood;
    }
    //calculate the probability at each site given read count(coverage)(negative-binomial distribution)
    public double logCoverageLikelihood(int c, double p, double r) {
        // negative binomial pmf
        double logCoverageLikelihood;
        logCoverageLikelihood = Gamma.logGamma(c + r) - Gamma.logGamma(r) - intLogGamma(c + 1) + r * Math.log(p) + c * Math.log(1 - p);
        return logCoverageLikelihood;
    }
    //get indices from propensities matrix by given genotype
    private static int[] getGenotypeIndices(int genotypeState) {
        int[][] indexTable = {
                {0,0},      //AA (AA and A_)
                {1,0,4},    //AC (AC, A_ and C_)
                {2,0,7},    //AG (AG, A_ and G_)
                {3,0,9},    //AT (AT, A_ and T_)
                {1,4,0},    //CA (CA, A_ and C_)
                {4,4},      //CC (CC and C_)
                {5,4,7},    //CG (CG, C_ and G_)
                {6,4,9},    //CT (CT, C_ and T_)
                {2,7,0},    //GA (GA, G_ and A_)
                {5,7,4},    //GC (GC, G_ and C_)
                {7,7},      //GG (GG and G_)
                {8,7,9},    //GT (GT, G_ and T_)
                {3,9,0},    //TA (TA, T_ and A_)
                {6,9,4},    //TC (TC, T_ and C_)
                {8,9,7},    //TG (TG, T_ and G_)
                {9,9},      //TT (TT and T_)
        };

        int[] indices = indexTable[genotypeState];
        return indices;
    }

    //Determining whether a genotype is homozygous or not
    private boolean homozygous(int genotype){
        return switch (genotype){
            case 0, 5, 10, 15 -> true;
            default -> false;
        };
    }

    //calculate the likelihood given read count (multinomial distribution)
    public double logLikelihoodDirichletMD(double w, int coverage, Double[] propensities, int[] readCountNumbers){
        double logLikelihood = logFFunction(coverage, w);
        for (int i = 0; i < readCountNumbers.length; i++) {
            logLikelihood = logLikelihood - logFFunction(readCountNumbers[i], w * propensities[i]);
        }
        return logLikelihood;
    }

    public double logFFunction(int coverage, double w){
        double result;
        if (coverage > 0){
            result = intLog(coverage) + Gamma.logGamma(w) + intLogGamma(coverage) - Gamma.logGamma(w + coverage);
            return result;
        } else return 0.0;
    }

    private double intLogGamma(int value){
        return logGammaCache.computeIfAbsent(value, Gamma::logGamma);
    }

    private double intLog(int value){
        return logCache.computeIfAbsent(value, Math::log);
    }


    @Override
    public void store() {
        super.store();
        /**
         * make a copy of current LogP's for each sequence
         * so that when the proposal is rejected it can be reversed
         **/
        System.arraycopy(currentLogPi, 0, storedLogPi, 0, storedLogPi.length);
    }

    @Override
    public void restore() {
        super.restore();

        /**
         * swap storedLogPi and currentLogPi, so that currentLogPi is now uptodate again
         */
        double[] tmp = storedLogPi;
        storedLogPi = currentLogPi;
        currentLogPi = tmp;
    }


}
