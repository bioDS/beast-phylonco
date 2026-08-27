package phylonco.beast.evolution.alignment;

import beast.base.core.BEASTInterface;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/** Round trip from the read-count NEXUS format lphy writes into a {@link ReadCountAlignment}. */
public class ReadCountNexusImporterTest {

    /** The format phylonco-lphy's ReadCountNexus writes, including its trailing ';' per row. */
    private static final String NEXUS = """
            #NEXUS

            begin taxa;
            \tdimensions ntax=2;
            \ttaxlabels cell1 cell2;
            end;

            begin characters;\tdimensions nchar=2;
            \tformat datatype=readCount;
            \tmatrix
            \t\tcell1 1:0:0:11,0:17:0:12;
            \t\tcell2 7:0:0:26,0:12:0:8;
            end;
            """;

    private File write(String name, String content) throws Exception {
        File f = File.createTempFile(name, ".nexus");
        f.deleteOnExit();
        Files.writeString(f.toPath(), content, StandardCharsets.UTF_8);
        return f;
    }

    @Test
    public void testImportsReadCountNexus() throws Exception {
        File file = write("readCounts", NEXUS);
        ReadCountNexusImporter importer = new ReadCountNexusImporter();

        assertTrue(importer.canHandleFile(file));
        List<BEASTInterface> loaded = importer.loadFile(file);
        assertEquals(1, loaded.size());

        ReadCountAlignment rc = (ReadCountAlignment) loaded.get(0);
        assertEquals(2, rc.getTaxaNumber());
        assertEquals(2, rc.getSiteNumber());
        assertEquals("cell1", rc.getTaxonName(0));
        assertEquals("cell2", rc.getTaxonName(1));
        assertArrayEquals(new int[]{1, 0, 0, 11}, rc.getReadCounts(0, 0));
        assertArrayEquals(new int[]{0, 17, 0, 12}, rc.getReadCounts(0, 1));
        assertArrayEquals(new int[]{7, 0, 0, 26}, rc.getReadCounts(1, 0));
        assertArrayEquals(new int[]{0, 12, 0, 8}, rc.getReadCounts(1, 1));
        assertEquals(33, rc.getMaxCoverage()); // cell2 site 0: 7 + 26
        // an Alignment BEAUti can hold as a partition
        assertEquals(2, rc.getPatternCount());
    }

    /** An ordinary genotype NEXUS must be left to the standard parser. */
    @Test
    public void testDoesNotClaimPlainNexus() throws Exception {
        URL url = getClass().getClassLoader().getResource("gt16ReadCountModel_A.nexus");
        assertTrue(url != null, "missing test resource");
        assertFalse(new ReadCountNexusImporter().canHandleFile(new File(url.toURI())));
    }

    @Test
    public void testRejectsFileWithoutReadCountDatatype() throws Exception {
        File file = write("notReadCounts", NEXUS.replace("datatype=readCount", "datatype=nucleotide"));
        assertThrows(IllegalArgumentException.class,
                () -> new ReadCountNexusImporter().loadFile(file));
    }
}
