package phylonco.beast.evolution.alignment;

import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.MCMC;
import beast.base.parser.XMLParser;
import beast.pkgmgmt.BEASTClassLoader;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import phylonco.beast.evolution.likelihood.ReadCountTreeLikelihood;

import java.io.File;
import java.net.URL;

import static beast.pkgmgmt.BEASTClassLoader.addServices;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * End-to-end XML checks for {@link ReadCountAlignment}: that BEAST parses it, that it produces the
 * same posterior as the legacy {@code ReadCount} + genotype-scaffold encoding of identical data, and
 * that an MCMC over it actually runs.
 */
public class ReadCountAlignmentXMLTest {

    private static final double DELTA = 1e-9;

    @BeforeAll
    public static void setUpClass() {
        BEASTClassLoader.initServices();
        addServices("version.xml");
    }

    private MCMC parse(String resource) throws Exception {
        URL url = getClass().getClassLoader().getResource(resource);
        assertTrue(url != null, "missing test resource " + resource);
        return (MCMC) new XMLParser().parseFile(new File(url.toURI()));
    }

    private static Distribution posteriorOf(MCMC mcmc) {
        return mcmc.posteriorInput.get();
    }

    /** MCMC sets up the State; a posterior evaluated before that returns NaN. */
    private static double initialLogP(MCMC mcmc) {
        return mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());
    }

    /** Walks the parsed posterior to the ReadCountTreeLikelihood's data alignment. */
    private static ReadCountAlignment dataOf(MCMC mcmc) {
        for (Distribution d : ((CompoundDistribution) mcmc.posteriorInput.get()).pDistributions.get()) {
            if (d instanceof ReadCountTreeLikelihood) {
                return (ReadCountAlignment) ((ReadCountTreeLikelihood) d).dataInput.get();
            }
        }
        throw new AssertionError("no ReadCountTreeLikelihood in the posterior");
    }

    /** The alignment parses from XML with the shape the rest of the model relies on. */
    @Test
    public void testParsesWithIdentityPatterns() throws Exception {
        ReadCountAlignment data = dataOf(parse("readCountAlignment.xml"));

        assertEquals(4, data.getTaxonCount());
        assertEquals(4, data.getTaxaNumber());
        assertEquals(6, data.getSiteCount());
        assertEquals(6, data.getSiteNumber());
        // identity patterns: pattern index == site index, no compression
        assertEquals(data.getSiteCount(), data.getPatternCount());
        for (int i = 0; i < data.getSiteCount(); i++) {
            assertEquals(i, data.getPatternIndex(i));
        }
        // the genotype datatype drives the tip-partial width
        assertEquals(10, data.getDataType().getStateCount());
        assertEquals("nucleotideDiploid10", data.getDataType().getTypeDescription());

        // read counts survive the round trip, in row order
        assertArrayEquals(new int[]{12, 0, 0, 0}, data.getReadCounts(0, 0));
        assertArrayEquals(new int[]{6, 5, 0, 0}, data.getReadCounts(0, 2));
        assertArrayEquals(new int[]{0, 0, 9, 0}, data.getReadCounts(3, 5));
        assertEquals("a", data.getTaxonName(0));
        assertEquals("d", data.getTaxonName(3));
        assertEquals(11, data.getCoverage(0, 2));
        assertEquals(13, data.getMaxCoverage()); // cell a, site 5: 7 + 6
    }

    /**
     * The new encoding and the legacy ReadCount + scaffold encoding of the same counts give the same
     * posterior, which is the real test that no data is lost or reordered by the Alignment wrapper.
     */
    @Test
    public void testSamePosteriorAsLegacyEncoding() throws Exception {
        double newLogP = initialLogP(parse("readCountAlignment.xml"));
        double legacyLogP = initialLogP(parse("readCountLegacy.xml"));

        assertTrue(Double.isFinite(newLogP), "posterior is not finite: " + newLogP);
        assertEquals(legacyLogP, newLogP, DELTA,
                "ReadCountAlignment posterior differs from the legacy ReadCount encoding");
    }

    /** A short MCMC over the read-count alignment runs to completion. */
    @Test
    public void testMCMCRuns() throws Exception {
        MCMC mcmc = parse("readCountAlignment.xml");
        File stateFile = File.createTempFile("readCountAlignmentTest", ".state");
        stateFile.deleteOnExit();
        mcmc.setStateFile(stateFile.getAbsolutePath(), false);
        mcmc.run();

        double logP = posteriorOf(mcmc).getCurrentLogP();
        assertTrue(Double.isFinite(logP), "posterior is not finite after the run: " + logP);
    }
}
