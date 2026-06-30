package phylonco.beast.evolution.readcountmodel;



import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import mutablealignment.MutableAlignment;

import org.apache.commons.numbers.gamma.LogGamma;
// replaces import org.apache.commons.math3.special.Gamma;
import phylonco.beast.evolution.datatype.NucleotideDiploid10;
import phylonco.beast.evolution.datatype.NucleotideDiploid16;
import phylonco.beast.evolution.datatype.ReadCount;

import java.util.List;
import java.util.Random;


public class LikelihoodReadCountModel extends Distribution {

    public Input<Alignment> alignmentInput = new Input<>("alignment",
            "optional genotype alignment; in the integrated model there is no genotype alignment and "
            + "the scaffold (taxa, site count, genotype datatype) is derived from the read counts "
            + "(datatype set via setDataType()) instead");
    public Input<ReadCount> readCountInput = new Input<>("readCount", "nucleotide read counts");

    // epsilon, allelic dropout, ... parameters
    public Input<RealScalarParam> epsilonInput = new Input<>("epsilon", "sequencing error");
    public Input<RealScalarParam> deltaInput = new Input<>("delta", "allelic dropout probability");
    public Input<RealScalarParam> tInput = new Input<>("t", "mean of allelic coverage");
    public Input<RealScalarParam> vInput = new Input<>("v", "variance of allelic coverage");
    public Input<RealVectorParam> sInput = new Input<>("s", "size factor of cell");
    public Input<RealScalarParam> w1Input = new Input<>("w1", "homozygous genotype overdispersion parameter of Dirichlet multinomial distribution");
    public Input<RealScalarParam> w2Input = new Input<>("w2", "heterogeneous genotype overdispersion parameter of Dirichlet multinomial distribution");

    // other parameters of read count model

    private RealScalarParam epsilon;
    private RealScalarParam delta;
    private RealScalarParam t;
    private RealScalarParam v;
    private RealVectorParam s;
    private RealScalarParam w1;
    private RealScalarParam w2;
    private Alignment alignment;
    private ReadCount readCount;
    private int numTaxa;   // number of cells (from alignment, or from readCount when integrated)
    private int numSites;  // number of sites (from alignment, or from readCount when integrated)
    private double[] negp1, negp2, negr1, negr2;
    private double[] wv;
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

    // Mapping from alignment taxon index to ReadCount taxon index.
    // These may differ when alignment sequences and ReadCount data rows
    // are in different orders (e.g. alignment order from XML vs alphabetical).
    private int[] alignToRCIndex;

    private int[][] gt16IndexTable;
    private int[][] gt10IndexTable;

    private boolean nbDirty = true;
    private boolean dirichletDirty = true;
    private boolean deltaDirty = true;

    DataType datatype;    //private double[] sv;

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

    private void initGt10IndexTable() {
        NucleotideDiploid10 datatype = new NucleotideDiploid10();
        // get genotype indices from data type NucleotideDiploid10
        gt10IndexTable = new int[10][];
        gt10IndexTable[0] = datatype.getIndices(new String[]{"AA", "A_"});
        gt10IndexTable[1] = datatype.getIndices(new String[]{"AC", "A_", "C_"});
        gt10IndexTable[2] = datatype.getIndices(new String[]{"AG", "A_", "G_"});
        gt10IndexTable[3] = datatype.getIndices(new String[]{"AT", "A_", "T_"});
        gt10IndexTable[4] = datatype.getIndices(new String[]{"CC", "C_"});
        gt10IndexTable[5] = datatype.getIndices(new String[]{"CG", "C_", "G_"});
        gt10IndexTable[6] = datatype.getIndices(new String[]{"CT", "C_", "T_"});
        gt10IndexTable[7] = datatype.getIndices(new String[]{"GG", "G_"});
        gt10IndexTable[8] = datatype.getIndices(new String[]{"GT", "G_", "T_"});
        gt10IndexTable[9] = datatype.getIndices(new String[]{"TT", "T_"});
    }

    private void initGt16IndexTable() {
        NucleotideDiploid16 datatype = new NucleotideDiploid16();
        // get genotype indices from data type NucleotideDiploid16
        gt16IndexTable = new int[22][];
        gt16IndexTable[0] = datatype.getIndices(new String[]{"AA", "A_"});
        gt16IndexTable[1] = datatype.getIndices(new String[]{"AC", "A_", "C_"});
        gt16IndexTable[2] = datatype.getIndices(new String[]{"AG", "A_", "G_"});
        gt16IndexTable[3] = datatype.getIndices(new String[]{"AT", "A_", "T_"});
        gt16IndexTable[4] = datatype.getIndices(new String[]{"CA", "C_", "A_"});
        gt16IndexTable[5] = datatype.getIndices(new String[]{"CC", "C_"});
        gt16IndexTable[6] = datatype.getIndices(new String[]{"CG", "C_", "G_"});
        gt16IndexTable[7] = datatype.getIndices(new String[]{"CT", "C_", "T_"});
        gt16IndexTable[8] = datatype.getIndices(new String[]{"GA", "G_", "A_"});
        gt16IndexTable[9] = datatype.getIndices(new String[]{"GC", "G_", "C_"});
        gt16IndexTable[10] = datatype.getIndices(new String[]{"GG", "G_"});
        gt16IndexTable[11] = datatype.getIndices(new String[]{"GT", "G_", "T_"});
        gt16IndexTable[12] = datatype.getIndices(new String[]{"TA", "T_", "A_"});
        gt16IndexTable[13] = datatype.getIndices(new String[]{"TC", "T_", "C_"});
        gt16IndexTable[14] = datatype.getIndices(new String[]{"TG", "T_", "G_"});
        gt16IndexTable[15] = datatype.getIndices(new String[]{"TT", "T_"});
        gt16IndexTable[16] = datatype.getIndices(new String[]{"AC", "A_", "C_"});  // AC or CA - M
        gt16IndexTable[17] = datatype.getIndices(new String[]{"AG", "A_", "G_"});  // AG or GA - R
        gt16IndexTable[18] = datatype.getIndices(new String[]{"AT", "A_", "T_"});  // AT or TA - W
        gt16IndexTable[19] = datatype.getIndices(new String[]{"CG", "C_", "G_"});  // CG or GC - S
        gt16IndexTable[20] = datatype.getIndices(new String[]{"CT", "C_", "T_"});  // CT or TC - Y
        gt16IndexTable[21] = datatype.getIndices(new String[]{"GT", "G_", "T_"});  // GT or TG - K
    }


    @Override
    public void initAndValidate() {
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

        String[] rcTaxaNames = readCount.getTaxaNames();
        if (alignment != null) {
            // standalone use: scaffold comes from the genotype alignment
            datatype = alignment.getDataType();
            buildIndexTable();
            numTaxa = alignment.getTaxonCount();
            numSites = alignment.getSiteCount();
            // Build mapping from alignment taxon index to ReadCount taxon index
            // (the alignment and ReadCount may have taxa in different orders).
            alignToRCIndex = new int[numTaxa];
            for (int i = 0; i < numTaxa; i++) {
                String alignTaxonName = alignment.getTaxaNames().get(i);
                alignToRCIndex[i] = i; // default: identity mapping
                for (int j = 0; j < rcTaxaNames.length; j++) {
                    if (alignTaxonName.equals(rcTaxaNames[j].trim())) {
                        alignToRCIndex[i] = j;
                        break;
                    }
                }
            }
        } else {
            // integrated model: there is no genotype alignment, so the scaffold (cells, sites,
            // taxon order) is the read-count data itself. The genotype datatype is supplied later
            // via setDataType() (from the GT10/GT16 substitution model in ReadCountTreeLikelihood).
            numTaxa = rcTaxaNames.length;
            numSites = readCount.getSiteNumber();
            alignToRCIndex = new int[numTaxa];
            for (int i = 0; i < numTaxa; i++) alignToRCIndex[i] = i; // identity
        }

        negp1 = new double[s.size()];
        negp2 = new double[s.size()];
        negr1 = new double[s.size()];
        negr2 = new double[s.size()];
        coverages = new int[numTaxa][numSites];
        for (int i = 0; i < numTaxa; i++) {
            int rcIdx = alignToRCIndex[i];
            for (int j = 0; j < numSites; j++) {
                for (int k = 0; k < 4; k++) {
                    coverages[i][j] += readCount.getReadCounts(rcIdx,j)[k];
                    if (readCount.getReadCounts(rcIdx,j)[k] > maxReadCount) {
                        maxReadCount = readCount.getReadCounts(rcIdx,j)[k];
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
            readDepthLogGamma[i] = LogGamma.value(i);
        }

        currentLogPi = new double[numTaxa];
        storedLogPi =  new double[numTaxa];
        rGammaLog = new double[2][s.size()];
        p1Log = new double[2][s.size()];
        p2Log = new double[2][s.size()];
        c_rLogGamma = new double[2][s.size()][maxReadDepth+1];
        c_wLogGamma = new double[2][maxReadDepth+1];
        rc_wPropLogGamma = new double[2][maxReadCount+1][10][4];
    }

    /** Builds the genotype-state -> nucleotide-pair propensity index table for the current datatype. */
    private void buildIndexTable() {
        if (datatype instanceof NucleotideDiploid10) {
            initGt10IndexTable();
        } else if (datatype instanceof NucleotideDiploid16) {
            initGt16IndexTable();
        }
    }

    /**
     * Sets the genotype datatype when there is no genotype alignment (integrated model). Called by
     * {@link phylonco.beast.evolution.likelihood.ReadCountTreeLikelihood} with the datatype implied
     * by the GT10/GT16 substitution model, and builds the genotype index table.
     */
    public void setDataType(DataType dataType) {
        this.datatype = dataType;
        buildIndexTable();
    }

    // calculate propensities matrix of dirichlet multinomial distribution(read count model)
    // and params of negative binomial distribution(coverage model)
    /**
     * Recompute cached values from current parameter values.
     * Public so GibbsSiteOperator can ensure caches are fresh before sampling
     * (they may be stale after a rejected parameter proposal), and so
     * ReadCountTreeLikelihood (a different package) can refresh the caches
     * before reading per-genotype likelihoods to build tip partials.
     */
    public void initialize() {
        double eps = epsilon.get();
        double del = delta.get();
        double tv = t.get();
        double vv = v.get();
        double[] sv = s.getValues();
        wv = new double[]{w1.get(), w2.get()};

        if (nbDirty) {
            updateNegativeBinomialCache(tv, vv, sv);
            nbDirty = false;
        }

        if (dirichletDirty) {
            updateDirichletCache(eps);
            dirichletDirty = false;
        }

        if (deltaDirty) {
            updateDeltaCache(del);
            deltaDirty = false;
        }
    }

    private void updateNegativeBinomialCache(double tv, double vv, double[] sv) {
        double mean1;
        double mean2;
        double variance1;
        double variance2;

        for (int i = 0; i < s.size(); i++) {
            mean1 = alpha[0] * tv * sv[i];
            mean2 = alpha[1] * tv * sv[i];

            variance1 = mean1 + Math.pow(alpha[0], 2) * vv * Math.pow(sv[i], 2);
            variance2 = mean2 + Math.pow(alpha[1], 2) * vv * Math.pow(sv[i], 2);

            negp1[i] = mean1 / variance1;
            negp2[i] = mean2 / variance2;
            negr1[i] = Math.pow(mean1, 2) / (variance1 - mean1);
            negr2[i] = Math.pow(mean2, 2) / (variance2 - mean2);

            rGammaLog[0][i] = LogGamma.value(negr1[i]);
            rGammaLog[1][i] = LogGamma.value(negr2[i]);

            p1Log[0][i] = Math.log(negp1[i]);
            p1Log[1][i] = Math.log(1 - negp1[i]);
            p2Log[0][i] = Math.log(negp2[i]);
            p2Log[1][i] = Math.log(1 - negp2[i]);

            for (int j = 0; j < maxReadDepth + 1; j++) {
                c_rLogGamma[0][i][j] = LogGamma.value(j + negr1[i]);
                c_rLogGamma[1][i][j] = LogGamma.value(j + negr2[i]);
            }
        }
    }

    private void updateDirichletCache(double eps) {
        double x0, x1, x2, x3;
        double[][] propensities;

        for (int i = 0; i < wPropensitiesLogGamma.length; i++) {
            x0 = LogGamma.value((1 - eps) * wv[i]);
            x1 = LogGamma.value((eps / 3) * wv[i]);
            x2 = LogGamma.value((0.5 - eps / 6) * wv[i]);
            x3 = LogGamma.value((eps / 6) * wv[i]);

            wPropensitiesLogGamma[i] = new double[][]{
                    {x0, x1, x1, x1},
                    {x2, x2, x3, x3},
                    {x2, x3, x2, x3},
                    {x2, x3, x3, x2},
                    {x1, x0, x1, x1},
                    {x3, x2, x2, x3},
                    {x3, x2, x3, x2},
                    {x1, x1, x0, x1},
                    {x3, x3, x2, x2},
                    {x1, x1, x1, x0},
            };
        }

        double y0, y1, y2, y3;
        for (int i = 0; i < rc_wPropLogGamma.length; i++) {
            y0 = (1 - eps) * wv[i];
            y1 = eps / 3 * wv[i];
            y2 = (0.5 - eps / 6) * wv[i];
            y3 = eps / 6 * wv[i];

            propensities = new double[][]{
                    {y0, y1, y1, y1},
                    {y2, y2, y3, y3},
                    {y2, y3, y2, y3},
                    {y2, y3, y3, y2},
                    {y1, y0, y1, y1},
                    {y3, y2, y2, y3},
                    {y3, y2, y3, y2},
                    {y1, y1, y0, y1},
                    {y3, y3, y2, y2},
                    {y1, y1, y1, y0},
            };

            for (int j = 0; j < rc_wPropLogGamma[i].length; j++) {
                for (int k = 0; k < rc_wPropLogGamma[i][j].length; k++) {
                    for (int l = 0; l < rc_wPropLogGamma[i][j][k].length; l++) {
                        rc_wPropLogGamma[i][j][k][l] =
                                LogGamma.value(propensities[k][l] + j);
                    }
                }
            }
        }

        for (int i = 0; i < c_wLogGamma.length; i++) {
            for (int j = 0; j < c_wLogGamma[i].length; j++) {
                c_wLogGamma[i][j] = LogGamma.value(wv[i] + j);
            }
        }

        wLogGamma[0] = LogGamma.value(wv[0]);
        wLogGamma[1] = LogGamma.value(wv[1]);
    }

    private void updateDeltaCache(double del) {
        deltaLog[0] = Math.log(del);
        deltaLog[1] = Math.log(1 - del);
    }

    /** Returns the mapping from alignment taxon index to ReadCount taxon index. */
    public int[] getAlignToRCIndex() {
        return alignToRCIndex;
    }

    /** Total read depth (summed over the 4 nucleotides) for an alignment taxon at a site. */
    public int getCoverage(int taxonIndex, int site) {
        return coverages[taxonIndex][site];
    }

    /** The ReadCount data backing this model. */
    public ReadCount getReadCount() {
        return readCount;
    }

    /** The genotype datatype (NucleotideDiploid16 or NucleotideDiploid10). */
    public DataType getDataType() {
        return datatype;
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
            int rcIdx = alignToRCIndex[i];
            for (int j = 0; j < alignment.getSiteCount(); j++) {///？
                // dirichlet multinomial pmf
                int patternIndex = alignment.getPatternIndex(j);
                int genotypeState = alignment.getPattern(i, patternIndex);
                int[] readCountNumbers = readCount.getReadCounts(rcIdx, j);
                logPi += logLiklihoodRC(genotypeState, readCountNumbers, coverages[i][j], rcIdx);
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
                int rcIdx = alignToRCIndex[i];
                for (int j = 0; j < mutableAlignment.getSiteCount(); j++) {///？
                    // dirichlet multinomial pmf
                    int patternIndex = mutableAlignment.getPatternIndex(j);
                    int genotypeState = mutableAlignment.getPattern(i, patternIndex);
                    int[] readCountNumbers = readCountInput.get().getReadCounts(rcIdx, j);
                    logPi += logLiklihoodRC(genotypeState, readCountNumbers, coverages[i][j], rcIdx);
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
                int rcIdx = alignToRCIndex[i];
                for (int j = 0; j < alignment.getSiteCount(); j++) {///？
                    // dirichlet multinomial pmf
                    int patternIndex = alignment.getPatternIndex(j);
                    int genotypeState = alignment.getPattern(i, patternIndex);
                    int[] readCountNumbers = readCount.getReadCounts(rcIdx, j);
                    logPi += logLiklihoodRC(genotypeState, readCountNumbers, coverages[i][j], rcIdx);
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
        int rcIdx = alignToRCIndex[taxonIndex];
        for (int j = 0; j < genotypeSequence.length; j++) {
            // dirichlet multinomial pmf
            int genotypeState = genotypeSequence[j];
            int[] readCountNumbers = readCount.getReadCounts(rcIdx, j);
            taxonLogP[j] = logLiklihoodRC(genotypeState, readCountNumbers, coverages[taxonIndex][j], rcIdx);
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
        double logPart0;
        double logPart1;
        double logPart2;
        double max;

        if (homozygous(genotypeState)) {
            logLikelihoodDirichletMDDiploid = logLikelihoodDirichletMD(0, coverage, readCountNumbers, wPropensitiesLogGamma[0][indices[0]], indices[0]);
            logCoverageLikelihoodDiploid = logCoverageLikelihood(coverage, negr2[taxonIndex], rGammaLog[1][taxonIndex], p2Log[0][taxonIndex], p2Log[1][taxonIndex], c_rLogGamma[1][taxonIndex][coverage]);
            logLikelihoodDirichletMDHaploid0 = logLikelihoodDirichletMD(0, coverage, readCountNumbers, wPropensitiesLogGamma[0][indices[1]], indices[1]);
            logCoverageLikelihoodHaploid = logCoverageLikelihood(coverage, negr1[taxonIndex], rGammaLog[0][taxonIndex], p1Log[0][taxonIndex], p1Log[1][taxonIndex], c_rLogGamma[0][taxonIndex][coverage]);
            logPart0 = logLikelihoodDirichletMDDiploid + logCoverageLikelihoodDiploid + deltaLog[1];
            logPart1 = logLikelihoodDirichletMDHaploid0 + logCoverageLikelihoodHaploid + deltaLog[0];
            max = Math.max(logPart0, logPart1);
            logLikelihood = max + Math.log(
                    Math.exp(logPart0 - max) + Math.exp(logPart1 - max)
            );
        } else {
            logLikelihoodDirichletMDDiploid = logLikelihoodDirichletMD(1, coverage, readCountNumbers, wPropensitiesLogGamma[1][indices[0]], indices[0]);
            logCoverageLikelihoodDiploid = logCoverageLikelihood(coverage, negr2[taxonIndex], rGammaLog[1][taxonIndex], p2Log[0][taxonIndex], p2Log[1][taxonIndex], c_rLogGamma[1][taxonIndex][coverage]);
            logLikelihoodDirichletMDHaploid0 = logLikelihoodDirichletMD(0, coverage, readCountNumbers, wPropensitiesLogGamma[0][indices[1]], indices[1]);
            logCoverageLikelihoodHaploid = logCoverageLikelihood(coverage, negr1[taxonIndex], rGammaLog[0][taxonIndex], p1Log[0][taxonIndex], p1Log[1][taxonIndex], c_rLogGamma[0][taxonIndex][coverage]);
            logLikelihoodDirichletMDHaploid1 = logLikelihoodDirichletMD(0, coverage, readCountNumbers, wPropensitiesLogGamma[0][indices[2]], indices[2]);
            logPart0 = logLikelihoodDirichletMDDiploid + logCoverageLikelihoodDiploid + deltaLog[1];
            logPart1 = log0_5 + logLikelihoodDirichletMDHaploid0 + logCoverageLikelihoodHaploid + deltaLog[0];
            logPart2 = log0_5 + logLikelihoodDirichletMDHaploid1 + logCoverageLikelihoodHaploid + deltaLog[0];
            max = Math.max(logPart0, Math.max(logPart1, logPart2));
            logLikelihood = max + Math.log(
                    Math.exp(logPart0 - max) + Math.exp(logPart1 - max) + Math.exp(logPart2 - max)
            );
        }
        return logLikelihood;
    }
    //calculate the probability at each site given read count(coverage)(negative-binomial distribution)
    public double logCoverageLikelihood(int c, double r, double rGammaLog, double pLog0, double pLog1, double c_rLogGamma) {
        // negative binomial pmf
        double logCoverageLikelihood;
        logCoverageLikelihood = c_rLogGamma - rGammaLog - readDepthLogGamma[c+1] + r * pLog0 + c * pLog1;
        return logCoverageLikelihood;
    }
    //get indices from propensities matrix by given genotype
    private int[] getGenotypeIndices(int genotypeState) {
        if (datatype instanceof NucleotideDiploid16) {
            return gt16IndexTable[genotypeState];
        } else if (datatype instanceof NucleotideDiploid10) {
            return gt10IndexTable[genotypeState];
        } else {
            throw new IllegalArgumentException("Unsupported genotype: " + datatype.getTypeDescription());
        }
    }

    //Determining whether a genotype is homozygous or not
    private boolean homozygous(int genotype){
        if (datatype instanceof NucleotideDiploid16) {
        return switch (genotype){
            case 0, 5, 10, 15 -> true;
            default -> false;
        };
        } else if (datatype instanceof NucleotideDiploid10) {
            return switch (genotype){
                case 0, 4, 7, 9 -> true;
                default -> false;
            };
        }else {
            throw new IllegalArgumentException("Unsupported genotype: " + datatype.getTypeDescription());
        }
    }

    //calculate the likelihood given read count (multinomial distribution)
    public double logLikelihoodDirichletMD( int wIndex, int coverage, int[] readCountNumbers, double[] wPropensitiesLogGamma, int index){
        double logLikelihood = 0.0;
        if (coverage > 0){
            logLikelihood = readDepthLog[coverage]
                    + wLogGamma[wIndex]
                    + readDepthLogGamma[coverage]
                    - c_wLogGamma[wIndex][coverage];
        }
        //double logLikelihood = logFFunctionCov(coverage, wIndex);

        for (int i = 0; i < readCountNumbers.length; i++) {
            int rc = readCountNumbers[i];
            if (rc > 0) {
                logLikelihood -=
                        readDepthLog[rc]
                                + wPropensitiesLogGamma[i]
                                + readDepthLogGamma[rc]
                                - rc_wPropLogGamma[wIndex][rc][index][i];
            }

            //logLikelihood -= logFFunctionRC(readCountNumbers[i], wPropensitiesLogGamma[i], rc_wPropLogGamma[wIndex][readCountNumbers[i]][index][i]);
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
    public boolean requiresRecalculation() {

        if (t.somethingIsDirty() ||
                v.somethingIsDirty() ||
                s.somethingIsDirty()) {
            nbDirty = true;
        }

        if (epsilon.somethingIsDirty() ||
                w1.somethingIsDirty() ||
                w2.somethingIsDirty()) {
            dirichletDirty = true;
        }

        if (delta.somethingIsDirty()) {
            deltaDirty = true;
        }

        return true;
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

        nbDirty = true;
        dirichletDirty = true;
        deltaDirty = true;
        initialize();

        double[] tmp = storedLogPi;
        storedLogPi = currentLogPi;
        currentLogPi = tmp;
    }


}
