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
    private double[][] propensities = new double[10][4];;
    private double[][][] wPropensitiesLogGamma = new double[2][10][4];
    private int[][] coverages;
    private final double[] alpha = new double[]{1.0, 2.0};
    private double [] currentLogPi, storedLogPi;
    private double[] wLogGamma = new double[2];
    private double[] deltaLog = new double[2];
    private double[][] p1Log;
    private double[][] p2Log;
    private double[][] rGammaLog;
    private final double log0_5 = Math.log(0.5);
    private int maxReadDepth = 0;
    private int maxReadCount = 0;
    private double[] readDepthLog;
    private double[] readDepthLogGamma;
    private double[][][] c_rLogGamma;
    private double[][][][] rc_wPropLogGamma;
    private double[][] c_wLogGamma;
    private static final int[][] indexTable = {
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
                    if (readCount.getReadCounts(i,j)[k] > maxReadCount) {
                        maxReadCount = readCount.getReadCounts(i,j)[k];
                    }
                }
                if (coverages[i][j] > maxReadDepth) {
                    maxReadDepth = coverages[i][j];
                }
            }
        }

        readDepthLog = new double[maxReadDepth+2];
        readDepthLogGamma = new double[maxReadDepth+2];
        for (int i = 1; i < maxReadDepth+2; i++) {
            readDepthLog[i] = Math.log(i);
            readDepthLogGamma[i] = Gamma.logGamma(i);
        }

        currentLogPi = new double[alignment.getTaxonCount()];
        storedLogPi =  new double[alignment.getTaxonCount()];
        rGammaLog = new double[2][s.getDimension()];
        p1Log = new double[2][s.getDimension()];
        p2Log = new double[2][s.getDimension()];
        c_rLogGamma = new double[2][s.getDimension()][maxReadDepth+1];
        c_wLogGamma = new double[2][maxReadDepth+1];
        rc_wPropLogGamma = new double[2][maxReadCount+1][10][4];
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
        w = new double[]{w1.getValue(), w2.getValue()};

        for (int i = 0; i < s.getDimension(); i++) {
            mean1 = alpha[0] * tv * sv[i];
            mean2 = alpha[1] * tv * sv[i];
            variance1 = mean1 + Math.pow(alpha[0], 2) * vv * Math.pow(sv[i], 2);
            variance2 = mean2 + Math.pow(alpha[1], 2) * vv * Math.pow(sv[i], 2);
            negp1[i] = mean1 / variance1;
            negp2[i] = mean2 / variance2;
            negr1[i] = Math.pow(mean1, 2) / (variance1 - mean1);
            negr2[i] = Math.pow(mean2, 2) / (variance2 - mean2);
            rGammaLog[0][i] = Gamma.logGamma(negr1[i]);
            rGammaLog[1][i] = Gamma.logGamma(negr2[i]);
            p1Log[0][i] = Math.log(negp1[i]);
            p1Log[1][i] = Math.log(1-negp1[i]);
            p2Log[0][i] = Math.log(negp2[i]);
            p2Log[1][i] = Math.log(1-negp2[i]);
            for (int j =1; j < maxReadDepth+1; j++) {
                c_rLogGamma[0][i][j] = Gamma.logGamma(j + negr1[i]);
                c_rLogGamma[1][i][j] = Gamma.logGamma(j + negr2[i]);
            }
        }

        double x0, x1, x2, x3;
        for (int i =0; i < wPropensitiesLogGamma.length; i++) {
            x0 = Gamma.logGamma((1 -eps)*w[i]);
            x1 = Gamma.logGamma((eps/3)*w[i]);
            x2 = Gamma.logGamma((0.5 - eps/6)*w[i]);
            x3 = Gamma.logGamma((eps/6)*w[i]);
            wPropensitiesLogGamma[i] = new double[][]{
                    {x0, x1, x1, x1},   // AA or A_ 0
                    {x2, x2, x3, x3},   // AC or CA 1
                    {x2, x3, x2, x3},   // AG or GA 2
                    {x2, x3, x3, x2},   // AT or TA 3
                    {x1, x0, x1, x1},   // CC or C_ 4
                    {x3, x2, x2, x3},   // CG or GC 5
                    {x3, x2, x3, x2},   // CT or TC 6
                    {x1, x1, x0, x1},   // GG or G_ 7
                    {x3, x3, x2, x2},   // GT or TG 8
                    {x1, x1, x1, x0},   // TT or T_ 9
            };
        }

        double y0, y1, y2, y3;
        for (int i = 0; i < rc_wPropLogGamma.length; i++) {
            y0 = (1 -eps)*w[i];
            y1 = eps/3*w[i];
            y2 = (0.5 - eps/6)*w[i];
            y3 = eps/6*w[i];
            propensities = new double[][]{
                    {y0, y1, y1, y1},   // AA or A_ 0
                    {y2, y2, y3, y3},   // AC or CA 1
                    {y2, y3, y2, y3},   // AG or GA 2
                    {y2, y3, y3, y2},   // AT or TA 3
                    {y1, y0, y1, y1},   // CC or C_ 4
                    {y3, y2, y2, y3},   // CG or GC 5
                    {y3, y2, y3, y2},   // CT or TC 6
                    {y1, y1, y0, y1},   // GG or G_ 7
                    {y3, y3, y2, y2},   // GT or TG 8
                    {y1, y1, y1, y0},   // TT or T_ 9
            };
            for (int j = 1; j < rc_wPropLogGamma[i].length; j++) {
                for (int k = 0; k < rc_wPropLogGamma[i][j].length; k++) {
                    for (int l = 0; l < rc_wPropLogGamma[i][j][k].length; l++) {
                        rc_wPropLogGamma[i][j][k][l] = Gamma.logGamma(propensities[k][l] + j);
                    }
                }
            }
        }

        for (int i = 0; i < c_wLogGamma.length; i++) {
            for (int j = 1; j < c_wLogGamma[i].length; j++) {
                c_wLogGamma[i][j] = Gamma.logGamma(w[i] + j);
            }
        }

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
            logLikelihoodDirichletMDDiploid = logLikelihoodDirichletMD(0, coverage, readCountNumbers, wPropensitiesLogGamma[0][indices[0]], indices[0]);
            logCoverageLikelihoodDiploid = logCoverageLikelihood(coverage,negp2[taxonIndex], negr2[taxonIndex], rGammaLog[1][taxonIndex], p2Log[0][taxonIndex], p2Log[1][taxonIndex], c_rLogGamma[1][taxonIndex][coverage]);
            logLikelihoodDirichletMDHaploid0 = logLikelihoodDirichletMD(0, coverage, readCountNumbers, wPropensitiesLogGamma[0][indices[1]], indices[1]);
            logCoverageLikelihoodHaploid = logCoverageLikelihood(coverage,negp1[taxonIndex], negr1[taxonIndex], rGammaLog[0][taxonIndex], p1Log[0][taxonIndex], p1Log[1][taxonIndex], c_rLogGamma[0][taxonIndex][coverage]);
            part0 = logLikelihoodDirichletMDDiploid + logCoverageLikelihoodDiploid + deltaLog[1];
            part1 = logLikelihoodDirichletMDHaploid0 + logCoverageLikelihoodHaploid + deltaLog[0];
            max = Math.max(part0, part1);
            if (part0 == max){
                logLikelihood = part0 + Math.log(1 + Math.exp(part1 - part0));
            }else {
                logLikelihood = part1 + Math.log(1 + Math.exp(part0 - part1));
            }
        } else {
            logLikelihoodDirichletMDDiploid = logLikelihoodDirichletMD(1, coverage, readCountNumbers, wPropensitiesLogGamma[1][indices[0]], indices[0]);
            logCoverageLikelihoodDiploid = logCoverageLikelihood(coverage, negp2[taxonIndex], negr2[taxonIndex], rGammaLog[1][taxonIndex], p2Log[0][taxonIndex], p2Log[1][taxonIndex], c_rLogGamma[1][taxonIndex][coverage]);
            logLikelihoodDirichletMDHaploid0 = logLikelihoodDirichletMD(0, coverage, readCountNumbers, wPropensitiesLogGamma[0][indices[1]], indices[1]);
            logCoverageLikelihoodHaploid = logCoverageLikelihood(coverage, negp1[taxonIndex], negr1[taxonIndex], rGammaLog[0][taxonIndex], p1Log[0][taxonIndex], p1Log[1][taxonIndex], c_rLogGamma[0][taxonIndex][coverage]);
            logLikelihoodDirichletMDHaploid1 = logLikelihoodDirichletMD(0, coverage, readCountNumbers, wPropensitiesLogGamma[0][indices[2]], indices[2]);
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
    public double logCoverageLikelihood(int c, double p, double r, double rGammaLog, double pLog0, double pLog1, double c_rLogGamma) {
        // negative binomial pmf
        double logCoverageLikelihood;
        logCoverageLikelihood = c_rLogGamma - rGammaLog - readDepthLogGamma[c+1] + r * pLog0 + c * pLog1;
        return logCoverageLikelihood;
    }
    //get indices from propensities matrix by given genotype
    private static int[] getGenotypeIndices(int genotypeState) {
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
    public double logLikelihoodDirichletMD( int wIndex, int coverage, int[] readCountNumbers, double[] wPropensitiesLogGamma, int index){
        double logLikelihood = logFFunctionCov(coverage, wIndex);
        for (int i = 0; i < readCountNumbers.length; i++) {
            logLikelihood = logLikelihood - logFFunctionRC(readCountNumbers[i], wPropensitiesLogGamma[i], rc_wPropLogGamma[wIndex][readCountNumbers[i]][index][i]);
        }
        return logLikelihood;
    }

    public double logFFunctionCov(int coverage, int wIndex){
        double result;
        if (coverage > 0){
            result = readDepthLog[coverage] + wLogGamma[wIndex] + readDepthLogGamma[coverage] - c_wLogGamma[wIndex][coverage];
            return result;
        } else return 0.0;
    }

    public double logFFunctionRC(int rc, double wPropensitiesLogGamma, double rc_wPropLogGamma){
        double result;
        if (rc > 0){
            result = readDepthLog[rc] + wPropensitiesLogGamma + readDepthLogGamma[rc] - rc_wPropLogGamma;
            return result;
        } else return 0.0;
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
