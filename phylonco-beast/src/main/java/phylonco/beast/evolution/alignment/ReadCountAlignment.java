package phylonco.beast.evolution.alignment;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.datatype.DataType;
import phylonco.beast.evolution.datatype.DiploidDataType;
import phylonco.beast.evolution.datatype.ReadCountMatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Single-cell nucleotide read counts as a BEAUti-visible {@link Alignment}: one row per cell, one
 * column per site, four counts (A,C,G,T) per cell/site.
 *
 * <p><b>Why an Alignment.</b> BEAUti can only make a partition out of an {@code Alignment}
 * ({@code BeautiDoc.alignments} is {@code List<Alignment>}; {@code BeautiAlignmentProvider} discards
 * anything else). Making the read-count data itself an Alignment also removes the identity-pattern
 * scaffold that {@code ReadCountTreeLikelihood} used to build internally and then drop from the XML.</p>
 *
 * <p><b>Datatype.</b> The declared datatype is the <i>genotype</i> datatype (nucleotideDiploid10 or
 * nucleotideDiploid16), because that is what the tree likelihood reads {@code getStateCount()} from
 * when sizing tip partials. The read counts live alongside, in {@link #readCounts}, and are reached
 * through {@link ReadCountMatrix}.</p>
 *
 * <p><b>The genotype states are cosmetic.</b> {@code counts} / {@code sitePatterns} hold a consensus
 * genotype call per cell/site, used for display, logging and the BEAUti alignment viewer only.
 * {@code TreeLikelihoodWithError} forces {@code useAmbiguities} and {@code useTipLikelihoods} true and
 * overrides {@code setPartials} to call {@code getLeafPartials}, so no tip state is ever read for the
 * likelihood — the partials come from the read counts.</p>
 *
 * <p><b>Patterns are identity.</b> {@link #calcPatterns(boolean)} is overridden so that pattern index
 * == site index, with no compression. {@code getLeafPartials} indexes read counts by pattern as if it
 * were a site, so column collapsing would silently misalign the data. (Setting
 * {@code usingTipLikelihoods} is not a substitute: the base method then makes one pattern per site but
 * still resolves {@code patternIndex} by binary search, so duplicate columns map to the wrong site.)</p>
 *
 * <p><b>StateNode.</b> {@code Alignment} extends {@code StateNode} but its store/restore/copy are
 * no-ops. This is constant data, so nothing is overridden and it must never appear in {@code <state>}.</p>
 */
@Description("Alignment of single-cell nucleotide read counts (A,C,G,T per cell per site)")
public class ReadCountAlignment extends Alignment implements ReadCountMatrix {

    /** Sequence value format: sites separated by ',', the four base counts by ':'. */
    public static final String SITE_SEPARATOR = ",";
    public static final String COUNT_SEPARATOR = ":";
    private static final String BASES = "ACGT";
    private static final int N_BASES = 4;

    final public Input<String> sitesIndexInput = new Input<>("sitesIndex",
            "optional space-separated genomic positions, one per site; reporting only");

    final public Input<Double> minorAlleleFreqInput = new Input<>("minorAlleleFreq",
            "minor-allele fraction at or above which the consensus genotype call is heterozygous. "
                    + "Affects display only — the likelihood reads the read counts, not the genotype states",
            0.2);

    // cells x sites x 4
    private int[][][] readCounts;
    // cells x sites, cached row sums
    private int[][] coverage;
    private int maxCoverage;
    private int maxReadCount;
    private String[] sitesIndex;
    private int missingState;

    public ReadCountAlignment() {
        // genotype datatype, not nucleotide: the tree likelihood sizes tip partials from it
        dataTypeInput.setValue("nucleotideDiploid16", this);
    }

    @Override
    public void initAndValidate() {
        rejectUnsupportedInputs();

        if (userDataTypeInput.get() != null) {
            m_dataType = userDataTypeInput.get();
        } else {
            initDataType(); // resolves dataTypeInput against the registered DataType services
        }
        if (!(m_dataType instanceof DiploidDataType)) {
            throw new IllegalArgumentException("ReadCountAlignment (" + getID() + ") requires a diploid "
                    + "genotype datatype (nucleotideDiploid10 or nucleotideDiploid16), got '"
                    + m_dataType.getTypeDescription() + "'");
        }
        missingState = m_dataType.stringToEncoding(String.valueOf(DataType.MISSING_CHAR)).get(0);

        sequences = sequenceInput.get();
        buildFromSequences(sequences, true);

        if (taxonSetInput.get() != null && taxonSetInput.get().getTaxonCount() > 0) {
            sortByTaxonSet(taxonSetInput.get());
        }
    }

    /**
     * Fills every field the {@link Alignment} contract requires: taxaNames, counts, stateCounts,
     * sitePatterns, patternWeight, patternIndex, maxStateCount. Re-runnable — BEAUti and
     * {@link #sortByTaxonSet} reorder the sequence list after construction.
     */
    private void buildFromSequences(List<Sequence> seqs, boolean log) {
        if (seqs == null || seqs.isEmpty()) {
            throw new IllegalArgumentException("ReadCountAlignment (" + getID() + ") has no sequences");
        }
        taxaNames.clear();
        counts.clear();
        stateCounts.clear();
        tipLikelihoods.clear();

        int nTaxa = seqs.size();
        int nSites = -1;
        readCounts = new int[nTaxa][][];
        for (int i = 0; i < nTaxa; i++) {
            Sequence seq = seqs.get(i);
            String taxon = seq.getTaxon().trim();
            if (taxaNames.contains(taxon)) {
                throw new IllegalArgumentException("Duplicate cell in ReadCountAlignment ("
                        + getID() + "): " + taxon);
            }
            int[][] row = parseCounts(seq.getData(), taxon);
            if (nSites < 0) {
                nSites = row.length;
            } else if (row.length != nSites) {
                throw new IllegalArgumentException("Cell " + taxon + " has " + row.length
                        + " sites, expected " + nSites);
            }
            readCounts[i] = row;
            taxaNames.add(taxon);
            stateCounts.add(m_dataType.getStateCount());
            tipLikelihoods.add(null);
        }

        coverage = new int[nTaxa][nSites];
        maxCoverage = 0;
        maxReadCount = 0;
        for (int i = 0; i < nTaxa; i++) {
            List<Integer> consensus = new ArrayList<>(nSites);
            for (int j = 0; j < nSites; j++) {
                int[] rc = readCounts[i][j];
                int total = 0;
                for (int k = 0; k < N_BASES; k++) {
                    total += rc[k];
                    maxReadCount = Math.max(maxReadCount, rc[k]);
                }
                coverage[i][j] = total;
                maxCoverage = Math.max(maxCoverage, total);
                consensus.add(consensusState(rc, total));
            }
            counts.add(consensus);
        }

        parseSitesIndex(nSites);
        isAscertained = false;
        calcPatterns(log);
    }

    private int[][] parseCounts(String value, String taxon) {
        String data = value.replaceAll("\\s", "");
        String[] sites = data.split(SITE_SEPARATOR);
        int[][] row = new int[sites.length][N_BASES];
        for (int j = 0; j < sites.length; j++) {
            String[] fields = sites[j].split(COUNT_SEPARATOR);
            if (fields.length != N_BASES) {
                throw new IllegalArgumentException("Cell " + taxon + " site " + j + ": expected "
                        + N_BASES + " counts (A:C:G:T), got '" + sites[j] + "'");
            }
            for (int k = 0; k < N_BASES; k++) {
                try {
                    row[j][k] = Integer.parseInt(fields[k]);
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException("Cell " + taxon + " site " + j
                            + ": non-integer read count '" + fields[k] + "'");
                }
                if (row[j][k] < 0) {
                    throw new IllegalArgumentException("Cell " + taxon + " site " + j
                            + ": negative read count " + row[j][k]);
                }
            }
        }
        return row;
    }

    private void parseSitesIndex(int nSites) {
        if (sitesIndexInput.get() == null) {
            sitesIndex = null;
            return;
        }
        sitesIndex = sitesIndexInput.get().trim().split("\\s+");
        if (sitesIndex.length != nSites) {
            throw new IllegalArgumentException("sitesIndex has " + sitesIndex.length
                    + " entries but the alignment has " + nSites + " sites");
        }
    }

    /**
     * Majority genotype call: the most abundant base, paired with the second-most abundant when that
     * one reaches {@code minorAlleleFreq} of the coverage. Zero coverage gives the missing state.
     * Alleles are emitted in canonical (sorted) order so a phased datatype gets a deterministic call.
     */
    private int consensusState(int[] rc, int total) {
        if (total == 0) {
            return missingState;
        }
        int first = 0;
        for (int k = 1; k < N_BASES; k++) {
            if (rc[k] > rc[first]) first = k;
        }
        int second = -1;
        for (int k = 0; k < N_BASES; k++) {
            if (k == first) continue;
            if (second < 0 || rc[k] > rc[second]) second = k;
        }
        char a = BASES.charAt(first);
        char b = rc[second] >= minorAlleleFreqInput.get() * total ? BASES.charAt(second) : a;
        String genotype = a <= b ? "" + a + b : "" + b + a;
        return ((DiploidDataType) m_dataType).getIndex(genotype);
    }

    /**
     * Identity patterns: one pattern per site, in site order, no compression. Mirrors
     * {@code MutableAlignment.calcPatterns}, for the reason given in the class javadoc.
     */
    @Override
    protected void calcPatterns(boolean log) {
        int taxonCount = counts.size();
        int siteCount = counts.get(0).size();

        sitePatterns = new int[siteCount][taxonCount];
        for (int i = 0; i < taxonCount; i++) {
            List<Integer> sites = counts.get(i);
            for (int j = 0; j < siteCount; j++) {
                sitePatterns[j][i] = sites.get(j);
            }
        }
        patternWeight = new int[siteCount];
        Arrays.fill(patternWeight, 1);
        patternIndex = new int[siteCount];
        for (int i = 0; i < siteCount; i++) {
            patternIndex[i] = i;
        }
        maxStateCount = m_dataType.getStateCount();

        if (log) {
            Log.info.println(getID() + ": " + taxonCount + " cells x " + siteCount + " sites, "
                    + "datatype " + m_dataType.getTypeDescription() + ", max coverage " + maxCoverage);
        }
    }

    /**
     * The base implementation rebuilds via the private {@code initializeWithSequenceList}, which would
     * try to parse the count strings through {@code Sequence.getSequence(dataType)}. Reorder and
     * rebuild here instead.
     */
    @Override
    public void sortByTaxonSet(TaxonSet toSortBy) {
        List<Sequence> sorted = new ArrayList<>(sequences);
        Collections.sort(sorted, (o1, o2) -> Integer.compare(
                toSortBy.getTaxonIndex(o1.getTaxon()), toSortBy.getTaxonIndex(o2.getTaxon())));
        sequences = sorted;
        buildFromSequences(sorted, false);
    }

    /**
     * Switches the genotype datatype (e.g. when BEAUti changes the substitution model from GT16 to
     * GT10) and recomputes the consensus calls and pattern block.
     */
    public void setGenotypeDataType(DataType dataType) {
        if (!(dataType instanceof DiploidDataType)) {
            throw new IllegalArgumentException("expected a diploid genotype datatype, got "
                    + dataType.getTypeDescription());
        }
        if (dataType.getStateCount() == m_dataType.getStateCount()) {
            return;
        }
        m_dataType = dataType;
        missingState = m_dataType.stringToEncoding(String.valueOf(DataType.MISSING_CHAR)).get(0);
        buildFromSequences(sequences, false);
    }

    private void rejectUnsupportedInputs() {
        if (Boolean.TRUE.equals(stripInvariantSitesInput.get())) {
            throw new IllegalArgumentException("ReadCountAlignment does not support strip=true "
                    + "(the read-count mapping requires one pattern per site)");
        }
        if (siteWeightsInput.get() != null) {
            throw new IllegalArgumentException("ReadCountAlignment does not support site weights");
        }
        if (Boolean.TRUE.equals(isAscertainedInput.get())) {
            throw new IllegalArgumentException("ReadCountAlignment does not support ascertainment correction");
        }
    }

    // ---- ReadCountMatrix ----

    @Override
    public int[] getReadCounts(int taxon, int site) {
        return readCounts[taxon][site];
    }

    @Override
    public int getTaxaNumber() {
        return getTaxonCount();
    }

    @Override
    public int getSiteNumber() {
        return getSiteCount();
    }

    @Override
    public String getTaxonName(int taxon) {
        return taxaNames.get(taxon);
    }

    // ---- read-count extras ----

    /** Total depth at one cell and site (cached row sum). */
    public int getCoverage(int taxon, int site) {
        return coverage[taxon][site];
    }

    /** Largest total depth over all cells and sites; sizes the model's log/logGamma lookup tables. */
    public int getMaxCoverage() {
        return maxCoverage;
    }

    /** Largest single-base count over all cells and sites. */
    public int getMaxReadCount() {
        return maxReadCount;
    }

    /** Genomic position of a site, or its index as a string when no sitesIndex was given. */
    public String getSiteIndex(int site) {
        return sitesIndex == null ? Integer.toString(site) : sitesIndex[site];
    }
}
