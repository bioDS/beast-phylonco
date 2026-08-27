package phylonco.beast.evolution.datatype;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;

/**
 * A cells x sites matrix of nucleotide read counts.
 *
 * <p>Deliberately minimal, and deliberately shaped so that {@link ReadCount} satisfies it without any
 * change to its methods — existing XML keeps working while new XML can use
 * {@code phylonco.beast.evolution.alignment.ReadCountAlignment}, which carries the same data as a
 * BEAUti-visible {@code Alignment}.</p>
 *
 * <p>Note there is deliberately no {@code getTaxaNames()} returning {@code String[]}: an
 * {@code Alignment} already declares {@code List<String> getTaxaNames()}, and the two erase to the
 * same signature, so no class could implement both. Use {@link #getTaxonName(int)} or the
 * {@link #getTaxonNames()} default.</p>
 */
@Description("A cells x sites matrix of nucleotide read counts")
public interface ReadCountMatrix extends BEASTInterface {

    /** @return the A,C,G,T counts at one cell and site; length 4, never null. */
    int[] getReadCounts(int taxon, int site);

    /** @return number of cells (rows). */
    int getTaxaNumber();

    /** @return number of sites (columns). */
    int getSiteNumber();

    /** @return name of the cell at the given row index. */
    String getTaxonName(int taxon);

    /** @return the cell names in row order. */
    default String[] getTaxonNames() {
        String[] names = new String[getTaxaNumber()];
        for (int i = 0; i < names.length; i++) {
            names[i] = getTaxonName(i);
        }
        return names;
    }
}
