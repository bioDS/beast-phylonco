package phylonco.beast.evolution.readcountmodel;



import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.inference.CalculationNode;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
import mutablealignment.MutableAlignment;

import org.apache.commons.numbers.gamma.LogGamma;
// replaces import org.apache.commons.math3.special.Gamma;
import phylonco.beast.evolution.datatype.NucleotideDiploid10;
import phylonco.beast.evolution.datatype.NucleotideDiploid16;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.beast.evolution.datatype.ReadCountMatrix;

import java.util.List;
import java.util.Random;


public class LikelihoodReadCountModel extends Distribution {

    public Input<Alignment> alignmentInput = new Input<>("alignment",
            "optional genotype alignment; in the integrated model there is no genotype alignment and "
            + "the scaffold (taxa, site count, genotype datatype) is derived from the read counts "
            + "(datatype set via setDataType()) instead");
    public Input<ReadCountMatrix> readCountInput = new Input<>("readCount", "nucleotide read counts");

    // epsilon, allelic dropout, ... parameters
    public Input<RealScalar> epsilonInput = new Input<>("epsilon", "sequencing error");
    public Input<RealScalar> deltaInput = new Input<>("delta", "allelic dropout probability");
    public Input<RealScalar> tInput = new Input<>("t", "mean of allelic coverage");
    public Input<RealScalar> vInput = new Input<>("v", "variance of allelic coverage");
    public Input<RealVector> sInput = new Input<>("s", "size factor of cell");
    public Input<RealScalar> w1Input = new Input<>("w1", "homozygous genotype overdispersion parameter of Dirichlet multinomial distribution");
    public Input<RealScalar> w2Input = new Input<>("w2", "heterogeneous genotype overdispersion parameter of Dirichlet multinomial distribution");

    // other parameters of read count model

    private RealScalar epsilon;
    private RealScalar delta;
    private RealScalar t;
    private RealScalar v;
    private RealVector s;
    private RealScalar w1;
    private RealScalar w2;
    private Alignment alignment;
    private ReadCountMatrix readCount;
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
        // wPropensitiesLogGamma / rc_wPropLogGamma are indexed by the fixed 10-row canonical,
        // phase-independent genotype-pair order (AA,AC,AG,AT,CC,CG,CT,GG,GT,TT — see
        // updateDirichletCache()), which matches NucleotideDiploid10's own state numbering.
        // GT16's own getIndex() additionally distinguishes phase (e.g. "GA" != "AG") and has
        // heterozygous ambiguity codes, so its raw state indices (0-21) cannot be used to index
        // that 10-row table directly (e.g. "G_" resolves to GT16's GG=10, overflowing it).
        // Every phased/ambiguous genotype below is therefore resolved through NucleotideDiploid10
        // using its canonical (alphabetically-ordered) pair string instead.
        NucleotideDiploid10 datatype = new NucleotideDiploid10();
        // get genotype indices from data type NucleotideDiploid16
        gt16IndexTable = new int[22][];
        gt16IndexTable[0] = datatype.getIndices(new String[]{"AA", "A_"});
        gt16IndexTable[1] = datatype.getIndices(new String[]{"AC", "A_", "C_"});
        gt16IndexTable[2] = datatype.getIndices(new String[]{"AG", "A_", "G_"});
        gt16IndexTable[3] = datatype.getIndices(new String[]{"AT", "A_", "T_"});
        gt16IndexTable[4] = datatype.getIndices(new String[]{"AC", "C_", "A_"});  // CA -> canonical AC
        gt16IndexTable[5] = datatype.getIndices(new String[]{"CC", "C_"});
        gt16IndexTable[6] = datatype.getIndices(new String[]{"CG", "C_", "G_"});
        gt16IndexTable[7] = datatype.getIndices(new String[]{"CT", "C_", "T_"});
        gt16IndexTable[8] = datatype.getIndices(new String[]{"AG", "G_", "A_"});  // GA -> canonical AG
        gt16IndexTable[9] = datatype.getIndices(new String[]{"CG", "G_", "C_"});  // GC -> canonical CG
        gt16IndexTable[10] = datatype.getIndices(new String[]{"GG", "G_"});
        gt16IndexTable[11] = datatype.getIndices(new String[]{"GT", "G_", "T_"});
        gt16IndexTable[12] = datatype.getIndices(new String[]{"AT", "T_", "A_"}); // TA -> canonical AT
        gt16IndexTable[13] = datatype.getIndices(new String[]{"CT", "T_", "C_"}); // TC -> canonical CT
        gt16IndexTable[14] = datatype.getIndices(new String[]{"GT", "T_", "G_"}); // TG -> canonical GT
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

        String[] rcTaxaNames = readCount.getTaxonNames();
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

        resolveSizeFactorDimension();

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
        covA = new double[s.size()][maxReadDepth + 1];
        covB = new double[s.size()][maxReadDepth + 1];
        covX = new double[s.size()][maxReadDepth + 1];
        covStamp = new int[s.size()][maxReadDepth + 1];
    }

    /** Builds the genotype-state -> nucleotide-pair propensity index table for the current datatype. */
    private void buildIndexTable() {
        termPairs = null;
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
    /**
     * {@code s} carries one size factor per cell. A BEAUti template cannot know the cell count when
     * it writes the parameter, so a dimension-1 {@code s} is broadcast to one entry per cell here.
     * Any other mismatch is an error: left alone it surfaces much later as an
     * ArrayIndexOutOfBoundsException inside the likelihood rather than as a usable message.
     */
    private void resolveSizeFactorDimension() {
        if (s.size() == numTaxa) {
            return;
        }
        if (s.size() == 1 && s instanceof RealVectorParam) {
            RealVectorParam sParam = (RealVectorParam) s;
            sParam.dimensionInput.setValue(numTaxa, sParam);
            sParam.initAndValidate();
            return;
        }
        throw new IllegalArgumentException("s (size factor) has dimension " + s.size()
                + " but there are " + numTaxa + " cells; give s one value per cell, or a single "
                + "value to be used for every cell");
    }

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
        double[] sv = ((RealVectorParam) s).getValues();
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

    /** Runs body(from, to) over a range [0, n), possibly split across threads (default: serially). */
    @FunctionalInterface
    public interface RangeRunner {
        void run(int n, java.util.function.BiConsumer<Integer, Integer> body);
    }

    private RangeRunner rangeRunner = (n, body) -> body.accept(0, n);

    /**
     * Lets the cache rebuilds split their loops (by cell for the coverage tables, by read count for
     * the Dirichlet tables) across threads; every entry is computed as before, so the values do not
     * change.
     */
    public void setRangeRunner(RangeRunner runner) {
        rangeRunner = runner;
    }

    private void updateNegativeBinomialCache(double tv, double vv, double[] sv) {
        coverageStamp++;
        double mean1;
        double mean2;
        double variance1;
        double variance2;

        rangeRunner.run(s.size(), (from, to) -> negativeBinomialCacheRange(tv, vv, sv, from, to));
    }

    private void negativeBinomialCacheRange(double tv, double vv, double[] sv, int from, int to) {
        double mean1, mean2, variance1, variance2;
        for (int i = from; i < to; i++) {
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
            y0 = (1 - eps) * wv[i];
            y1 = eps / 3 * wv[i];
            y2 = (0.5 - eps / 6) * wv[i];
            y3 = eps / 6 * wv[i];

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

            final int wi = i;
            final double[][] prop = propensities;
            rangeRunner.run(rc_wPropLogGamma[wi].length, (from, to) -> {
                for (int j = from; j < to; j++) {
                    for (int k = 0; k < rc_wPropLogGamma[wi][j].length; k++) {
                        for (int l = 0; l < rc_wPropLogGamma[wi][j][k].length; l++) {
                            rc_wPropLogGamma[wi][j][k][l] =
                                    LogGamma.value(prop[k][l] + j);
                        }
                    }
                }
            });
        }

        for (int i = 0; i < c_wLogGamma.length; i++) {
            final int wi = i;
            rangeRunner.run(c_wLogGamma[wi].length, (from, to) -> {
                for (int j = from; j < to; j++) {
                    c_wLogGamma[wi][j] = LogGamma.value(wv[wi] + j);
                }
            });
        }

        wLogGamma[0] = LogGamma.value(wv[0]);
        wLogGamma[1] = LogGamma.value(wv[1]);
    }

    private void updateDeltaCache(double del) {
        coverageStamp++;
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
    public ReadCountMatrix getReadCount() {
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
    // ------------------------------------------------------------------ GT10 tip likelihoods by terms
    //
    // logLiklihoodRC evaluates, for each genotype, two or three Dirichlet-multinomial terms (the
    // diploid one and the haploid ones under allelic dropout). Over the ten GT10 genotypes these are
    // only a few distinct (overdispersion, propensity row) pairs, and they depend on epsilon, w1 and
    // w2 alone, not on the coverage parameters t, v, s or on delta. ReadCountTreeLikelihood caches
    // them per cell and site (dirichletTerms) and rebuilds the ten likelihoods from the cache
    // (logLikelihoodsFromTerms) with the same arithmetic as logLiklihoodRC, so the results are
    // identical.

    private int[][] termPairs;    // distinct {wIndex, propensity row}
    private int[][] genotypeTerms; // genotype -> indices into termPairs: diploid, haploid(s)

    private void ensureTermPairs() {
        if (termPairs != null) return;
        if (!(datatype instanceof NucleotideDiploid10)) {
            throw new IllegalStateException("tip likelihoods by terms are implemented for GT10 only");
        }
        java.util.List<int[]> pairs = new java.util.ArrayList<>();
        genotypeTerms = new int[10][];
        for (int g = 0; g < 10; g++) {
            int[] idx = gt10IndexTable[g];
            boolean hom = homozygous(g);
            int n = hom ? 2 : 3;
            genotypeTerms[g] = new int[n];
            for (int k = 0; k < n; k++) {
                int[] pair = {k == 0 && !hom ? 1 : 0, idx[k]};
                int found = -1;
                for (int q = 0; q < pairs.size(); q++) {
                    if (pairs.get(q)[0] == pair[0] && pairs.get(q)[1] == pair[1]) found = q;
                }
                if (found < 0) { pairs.add(pair); found = pairs.size() - 1; }
                genotypeTerms[g][k] = found;
            }
        }
        termPairs = pairs.toArray(new int[0][]);
    }

    /** Number of distinct Dirichlet-multinomial terms per cell and site (GT10). */
    public int dirichletTermCount() {
        ensureTermPairs();
        return termPairs.length;
    }

    /**
     * The distinct Dirichlet-multinomial log terms of one cell and site, written to out[off ...].
     * They depend on epsilon, w1 and w2 only; call initialize() first.
     */
    public void dirichletTerms(int[] readCountNumbers, int coverage, double[] out, int off) {
        ensureTermPairs();
        for (int q = 0; q < termPairs.length; q++) {
            int w = termPairs[q][0], row = termPairs[q][1];
            out[off + q] = logLikelihoodDirichletMD(w, coverage, readCountNumbers, wPropensitiesLogGamma[w][row], row);
        }
    }

    /**
     * logLiklihoodRC for each of the ten GT10 genotypes, from the cell and site's Dirichlet terms
     * (dirichletTerms at the current epsilon, w1, w2) and the current coverage and dropout
     * parameters; the same arithmetic as logLiklihoodRC, so the values are identical.
     */
    public void logLikelihoodsFromTerms(double[] terms, int off, int coverage, int taxonIndex, double[] out) {
        ensureTermPairs();
        double covDiploid = logCoverageLikelihood(coverage, negr2[taxonIndex], rGammaLog[1][taxonIndex],
                p2Log[0][taxonIndex], p2Log[1][taxonIndex], c_rLogGamma[1][taxonIndex][coverage]);
        double covHaploid = logCoverageLikelihood(coverage, negr1[taxonIndex], rGammaLog[0][taxonIndex],
                p1Log[0][taxonIndex], p1Log[1][taxonIndex], c_rLogGamma[0][taxonIndex][coverage]);
        for (int g = 0; g < 10; g++) {
            int[] k = genotypeTerms[g];
            double logPart0 = terms[off + k[0]] + covDiploid + deltaLog[1];
            if (k.length == 2) {
                double logPart1 = terms[off + k[1]] + covHaploid + deltaLog[0];
                double max = Math.max(logPart0, logPart1);
                out[g] = max + Math.log(Math.exp(logPart0 - max) + Math.exp(logPart1 - max));
            } else {
                double logPart1 = log0_5 + terms[off + k[1]] + covHaploid + deltaLog[0];
                double logPart2 = log0_5 + terms[off + k[2]] + covHaploid + deltaLog[0];
                double max = Math.max(logPart0, Math.max(logPart1, logPart2));
                out[g] = max + Math.log(Math.exp(logPart0 - max) + Math.exp(logPart1 - max) + Math.exp(logPart2 - max));
            }
        }
    }

    /**
     * Linear-space form of the cached terms for tipPartialsFromTerms: out[off + q] = exp(terms[q] - max)
     * for the cell and site's terms (dirichletTerms), and out[off + count] = max.
     */
    public void linearTerms(double[] terms, int termOff, double[] out, int off) {
        int n = dirichletTermCount();
        double max = Double.NEGATIVE_INFINITY;
        for (int q = 0; q < n; q++) max = Math.max(max, terms[termOff + q]);
        for (int q = 0; q < n; q++) out[off + q] = Math.exp(terms[termOff + q] - max);
        out[off + n] = max;
    }

    // coverage and dropout weights by cell and depth: for a diploid read set exp(covDiploid +
    // log(1 - delta) - X), for a haploid one exp(covHaploid + log(delta) - X), and X; rebuilt lazily
    // after any change to t, v, s or delta
    private double[][] covA, covB, covX;     // [taxon][depth]
    private int[][] covStamp;                // entry valid when equal to coverageStamp
    private int coverageStamp = 1;

    // allocated in initAndValidate, so that threads filling different cells' rows share no state
    private void ensureCoverageWeight(int taxonIndex, int coverage) {
        if (covStamp[taxonIndex][coverage] != coverageStamp) {
            double lA = logCoverageLikelihood(coverage, negr2[taxonIndex], rGammaLog[1][taxonIndex],
                    p2Log[0][taxonIndex], p2Log[1][taxonIndex], c_rLogGamma[1][taxonIndex][coverage]) + deltaLog[1];
            double lB = logCoverageLikelihood(coverage, negr1[taxonIndex], rGammaLog[0][taxonIndex],
                    p1Log[0][taxonIndex], p1Log[1][taxonIndex], c_rLogGamma[0][taxonIndex][coverage]) + deltaLog[0];
            double x = Math.max(lA, lB);
            covA[taxonIndex][coverage] = Math.exp(lA - x);
            covB[taxonIndex][coverage] = Math.exp(lB - x);
            covX[taxonIndex][coverage] = x;
            covStamp[taxonIndex][coverage] = coverageStamp;
        }
    }

    /**
     * Max-normalised tip partials of one cell and site for the ten GT10 genotypes, written to
     * out[outOff ...], from its cached terms in log form (terms, dirichletTerms) and linear form
     * (linear, linearTerms); returns the log of the normalising factor, so that
     * log P(reads | g) = log(out[g]) + return value, equal to logLiklihoodRC up to rounding. A
     * cell and site whose terms would underflow in linear form is computed by logLikelihoodsFromTerms.
     */
    public double tipPartialsFromTerms(double[] terms, int termOff, double[] linear, int linOff, int coverage,
                                       int taxonIndex, double[] out, int outOff, double[] scratch) {
        ensureTermPairs();
        ensureCoverageWeight(taxonIndex, coverage);
        double a = covA[taxonIndex][coverage], b = covB[taxonIndex][coverage];
        double m = 0.0;
        for (int g = 0; g < 10; g++) {
            int[] k = genotypeTerms[g];
            double v = k.length == 2
                    ? a * linear[linOff + k[0]] + b * linear[linOff + k[1]]
                    : a * linear[linOff + k[0]] + 0.5 * b * (linear[linOff + k[1]] + linear[linOff + k[2]]);
            out[outOff + g] = v;
            if (v > m) m = v;
        }
        if (m > 1e-250) {
            double inv = 1.0 / m;
            for (int g = 0; g < 10; g++) out[outOff + g] *= inv;
            return Math.log(m) + linear[linOff + termPairs.length] + covX[taxonIndex][coverage];
        }
        // underflow in linear form: the log-space computation
        logLikelihoodsFromTerms(terms, termOff, coverage, taxonIndex, scratch);
        double max = Double.NEGATIVE_INFINITY;
        for (int g = 0; g < 10; g++) if (scratch[g] > max) max = scratch[g];
        for (int g = 0; g < 10; g++) out[outOff + g] = Math.exp(scratch[g] - max);
        return max;
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

        if (((CalculationNode) t).somethingIsDirty() ||
                ((CalculationNode) v).somethingIsDirty() ||
                ((CalculationNode) s).somethingIsDirty()) {
            nbDirty = true;
        }

        if (((CalculationNode) epsilon).somethingIsDirty() ||
                ((CalculationNode) w1).somethingIsDirty() ||
                ((CalculationNode) w2).somethingIsDirty()) {
            dirichletDirty = true;
        }

        if (((CalculationNode) delta).somethingIsDirty()) {
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
