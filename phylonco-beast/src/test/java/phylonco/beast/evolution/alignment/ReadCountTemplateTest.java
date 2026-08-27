package phylonco.beast.evolution.alignment;

import beast.base.core.BEASTInterface;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.MCMC;
import beast.base.parser.XMLParser;
import beast.pkgmgmt.BEASTClassLoader;
import beastfx.app.inputeditor.BeautiDoc;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;
import javafx.application.Platform;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import static beast.pkgmgmt.BEASTClassLoader.addServices;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Drives BEAUti headlessly over the read-count template: load the template, add a read-count
 * partition the way the import menu would, write the XML, then parse and run it in BEAST. This is
 * the check that the template's connectors actually produce a runnable analysis.
 */
public class ReadCountTemplateTest {

    private static final String NEXUS = """
            #NEXUS

            begin taxa;
            \tdimensions ntax=4;
            \ttaxlabels cell1 cell2 cell3 cell4;
            end;

            begin characters;\tdimensions nchar=6;
            \tformat datatype=readCount;
            \tmatrix
            \t\tcell1 12:0:0:0,0:10:0:0,6:5:0:0,0:0:9:0,0:0:0:11,7:0:6:0;
            \t\tcell2 10:0:0:0,0:9:0:0,8:0:0:0,0:0:7:0,0:0:0:9,0:0:12:0;
            \t\tcell3 9:0:0:0,0:11:0:0,5:6:0:0,0:0:8:0,0:0:0:10,6:0:5:0;
            \t\tcell4 11:0:0:0,0:8:0:0,7:0:0:0,0:0:10:0,0:0:0:8,0:0:9:0;
            end;
            """;

    @BeforeAll
    public static void setUpClass() {
        // BeautiDoc lives in beastfx, so the JavaFX toolkit has to exist even though this test
        // never shows a window
        System.setProperty("prism.order", "sw");
        try {
            Platform.startup(() -> { });
        } catch (IllegalStateException alreadyRunning) {
            // fine: the toolkit is up
        }
        BEASTClassLoader.initServices();
        addServices("version.xml");
    }

    private static File resourceFile(String name) throws Exception {
        URL url = ReadCountTemplateTest.class.getClassLoader().getResource(name);
        assertTrue(url != null, "missing resource " + name);
        return new File(url.toURI());
    }

    @Test
    public void testTemplateGeneratesRunnableXML() throws Exception {
        File nexus = File.createTempFile("readCounts", ".nexus");
        nexus.deleteOnExit();
        Files.writeString(nexus.toPath(), NEXUS, StandardCharsets.UTF_8);

        BeautiDoc doc = new BeautiDoc();
        String templateXML = doc.processTemplate(
                resourceFile("fxtemplates/ReadCountModel.xml").getAbsolutePath());
        doc.initialize(BeautiDoc.ActionOnExit.UNKNOWN, null, templateXML, "readCountTemplateTest");

        BEASTInterface imported = new ReadCountNexusImporter().loadFile(nexus).get(0);
        doc.addAlignmentWithSubnet((ReadCountAlignment) imported,
                doc.beautiConfig.partitionTemplate.get());

        File out = File.createTempFile("readCountBeauti", ".xml");
        out.deleteOnExit();
        doc.save(out);

        String xml = Files.readString(out.toPath(), StandardCharsets.UTF_8);
        assertTrue(xml.contains("ReadCountTreeLikelihood"), "no read-count tree likelihood in:\n" + xml);
        assertTrue(xml.contains("LikelihoodReadCountModel"), "no read-count model in:\n" + xml);
        assertTrue(xml.contains("ReadCountAlignment"), "read counts not written as data in:\n" + xml);

        // the tree prior and its operators are connected, so the tree is actually sampled
        assertTrue(xml.contains("YuleModel.t:"), "no tree prior in the XML");
        for (String op : new String[]{"YuleModelTreeScaler", "YuleModelUniformOperator",
                "YuleModelSubtreeSlide", "YuleModelNarrow", "YuleModelWide",
                "YuleModelWilsonBalding", "YuleBirthRateScaler"}) {
            assertTrue(xml.contains(op + ".t:"), "missing tree operator " + op);
        }

        // every read-count parameter reached the state, and got an operator
        for (String param : new String[]{"rcEpsilon", "rcDelta", "rcT", "rcV", "rcS", "rcW1", "rcW2"}) {
            assertTrue(xml.contains(param + ".s:"), "parameter " + param + " missing from the XML");
            assertTrue(xml.contains(param + "Scaler.s:"), "no operator for " + param);
        }

        MCMC mcmc = (MCMC) new XMLParser().parseFile(out);

        // the read-count model must be wired in only through the tree likelihood. If a connector
        // also put it in the posterior it would be counted twice.
        assertEquals(0, countReadCountModelsIn(mcmc.posteriorInput.get()),
                "LikelihoodReadCountModel is in the posterior as well as in the tree likelihood, "
                        + "so its contribution would be double counted");

        mcmc.chainLengthInput.setValue(200L, mcmc);
        mcmc.initAndValidate();
        File stateFile = File.createTempFile("readCountTemplate", ".state");
        stateFile.deleteOnExit();
        mcmc.setStateFile(stateFile.getAbsolutePath(), false);
        mcmc.run();

        double logP = mcmc.posteriorInput.get().getCurrentLogP();
        assertTrue(Double.isFinite(logP), "BEAUti-generated XML has a non-finite posterior: " + logP);
    }

    /** Counts LikelihoodReadCountModel instances reachable as posterior sub-distributions. */
    private static int countReadCountModelsIn(Distribution d) {
        if (d instanceof LikelihoodReadCountModel) {
            return 1;
        }
        int n = 0;
        if (d instanceof CompoundDistribution) {
            for (Distribution child : ((CompoundDistribution) d).pDistributions.get()) {
                n += countReadCountModelsIn(child);
            }
        }
        return n;
    }
}
