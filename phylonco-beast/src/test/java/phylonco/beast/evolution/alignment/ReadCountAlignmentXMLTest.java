package phylonco.beast.evolution.alignment;

import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.MCMC;
import beast.base.parser.XMLParser;
import beast.pkgmgmt.BEASTClassLoader;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import phylonco.beast.evolution.likelihood.ReadCountTreeLikelihood;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;

import java.io.File;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import static beast.pkgmgmt.BEASTClassLoader.addServices;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
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

    /** The base XML's four-per-cell size factor, as written into a temp copy of the XML. */
    private static final String S_FOUR =
            "<parameter id=\"s\" dimension=\"4\" spec=\"beast.base.spec.inference.parameter.RealVectorParam\" domain=\"PositiveReal\">\n"
            + "        1.04 1.04 1.02 1.01\n    </parameter>";

    /** Writes a copy of the base XML with the s parameter replaced, and parses it. */
    private MCMC parseWithSizeFactor(String replacement) throws Exception {
        URL url = getClass().getClassLoader().getResource("readCountAlignment.xml");
        String xml = Files.readString(new File(url.toURI()).toPath(), StandardCharsets.UTF_8);
        assertTrue(xml.contains(S_FOUR), "base XML no longer contains the expected s parameter");
        File out = File.createTempFile("readCountAlignmentS", ".xml");
        out.deleteOnExit();
        Files.writeString(out.toPath(), xml.replace(S_FOUR, replacement), StandardCharsets.UTF_8);
        return (MCMC) new XMLParser().parseFile(out);
    }

    private static LikelihoodReadCountModel readCountModelOf(MCMC mcmc) {
        for (Distribution d : ((CompoundDistribution) mcmc.posteriorInput.get()).pDistributions.get()) {
            if (d instanceof ReadCountTreeLikelihood) {
                return ((ReadCountTreeLikelihood) d).readCountModelInput.get();
            }
        }
        throw new AssertionError("no ReadCountTreeLikelihood in the posterior");
    }

    /**
     * A BEAUti template cannot know the cell count, so it writes a single size factor. That must be
     * broadcast to one entry per cell rather than silently mis-sizing the model's arrays.
     */
    @Test
    public void testSingleSizeFactorIsBroadcastToCells() throws Exception {
        MCMC mcmc = parseWithSizeFactor(
                "<parameter id=\"s\" spec=\"beast.base.spec.inference.parameter.RealVectorParam\" domain=\"PositiveReal\">1.04</parameter>");

        assertEquals(4, readCountModelOf(mcmc).sInput.get().size());
        assertTrue(Double.isFinite(initialLogP(mcmc)));
    }

    /** Any other dimension mismatch is a usable error, not an AIOOBE inside the likelihood. */
    @Test
    public void testWrongSizeFactorDimensionIsRejected() throws Exception {
        Exception e = assertThrows(Exception.class, () -> parseWithSizeFactor(
                "<parameter id=\"s\" dimension=\"3\" spec=\"beast.base.spec.inference.parameter.RealVectorParam\" domain=\"PositiveReal\">1.0 1.0 1.0</parameter>"));

        String message = "";
        for (Throwable t = e; t != null; t = t.getCause()) {
            message += t.getMessage() + " ";
        }
        assertTrue(message.contains("one value per cell"),
                "expected a size-factor dimension message, got: " + message);
    }

    /**
     * The substitution model is authoritative for the genotype state count. BEAUti lets the user
     * switch GT16 to GT10 in the Site Model tab after the partition exists, so a ReadCountAlignment
     * declaring the other datatype is reconciled rather than rejected.
     */
    @Test
    public void testSubstModelDrivesTheGenotypeDataType() throws Exception {
        URL url = getClass().getClassLoader().getResource("readCountAlignment.xml");
        String xml = Files.readString(new File(url.toURI()).toPath(), StandardCharsets.UTF_8);
        // declare 16 genotype states on the data while the substitution model stays GT10
        xml = xml.replace("nucleotideDiploid10", "nucleotideDiploid16");
        File out = File.createTempFile("readCountAlignmentDataType", ".xml");
        out.deleteOnExit();
        Files.writeString(out.toPath(), xml, StandardCharsets.UTF_8);

        MCMC mcmc = (MCMC) new XMLParser().parseFile(out);
        ReadCountAlignment data = dataOf(mcmc);

        assertEquals(10, data.getDataType().getStateCount(),
                "the GT10 substitution model should have driven the alignment back to 10 states");
        assertTrue(Double.isFinite(initialLogP(mcmc)));
    }
}
