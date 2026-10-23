package phylonco.beast.evolution.alignment;

import beast.base.core.BEASTInterface;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.MCMC;
import beast.base.parser.XMLParser;
import beast.pkgmgmt.BEASTClassLoader;
import beastfx.app.inputeditor.BeautiDoc;
import phylonco.beast.evolution.substitutionmodel.GT16;
import phylonco.beast.evolution.likelihood.ReadCountTreeLikelihood;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.evolution.substitutionmodel.Frequencies;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.domain.PositiveReal;
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
import static org.junit.jupiter.api.Assertions.assertFalse;
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

    private static File writeNexus() throws Exception {
        File nexus = File.createTempFile("readCounts", ".nexus");
        nexus.deleteOnExit();
        Files.writeString(nexus.toPath(), NEXUS, StandardCharsets.UTF_8);
        return nexus;
    }

    /** Loads the template and adds the read-count partition, as the import menu would. */
    private static BeautiDoc beautiDocFor(File nexus) throws Exception {
        BeautiDoc doc = new BeautiDoc();
        String templateXML = doc.processTemplate(
                resourceFile("fxtemplates/ReadCountModel.xml").getAbsolutePath());
        doc.initialize(BeautiDoc.ActionOnExit.UNKNOWN, null, templateXML, "readCountTemplateTest");

        BEASTInterface imported = new ReadCountNexusImporter().loadFile(nexus).get(0);
        doc.addAlignmentWithSubnet((ReadCountAlignment) imported,
                doc.beautiConfig.partitionTemplate.get());
        return doc;
    }

    @Test
    public void testTemplateGeneratesRunnableXML() throws Exception {
        File nexus = writeNexus();
        BeautiDoc doc = beautiDocFor(nexus);

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

        // every read-count parameter is estimated: a state node, with a prior and an operator
        String state = xml.substring(xml.indexOf("<state "), xml.indexOf("</state>"));
        for (String param : new String[]{"rcEpsilon", "rcDelta", "rcT", "rcV", "rcS", "rcW1", "rcW2"}) {
            assertTrue(state.contains("id=\"" + param + ".s:"),
                    param + " is not a state node, so it would be held fixed:\n" + state);
            assertTrue(xml.contains(param + "Prior.s:"), "no prior for " + param);
            assertTrue(xml.contains(param + "Scaler.s:"), "no operator for " + param);
        }
        // the single size factor in the template was broadcast to one per cell
        assertTrue(state.contains("id=\"rcS.s:") && state.contains("dimension=\"4\""),
                "rcS was not resized to the cell count:\n" + state);

        // the alignment is constant data: it must not be registered as a state node
        assertFalse(state.contains("<stateNode id=\"" + doc.alignments.get(0).getID()),
                "the read-count alignment was added to the state");

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

    /**
     * The analysis BEAUti generates must be the same model a hand-built one is. This rebuilds the
     * read-count tree likelihood independently, from the separately imported alignment and the
     * template's documented parameter values, over the tree BEAUti produced, and checks the two
     * agree. A mis-wired template (wrong data on the likelihood, a missing readCountModel, the
     * wrong site model) would show up here as a different likelihood, where a finiteness check
     * alone would pass.
     */
    @Test
    public void testGeneratedLikelihoodMatchesAHandBuiltOne() throws Exception {
        File nexus = writeNexus();
        BeautiDoc doc = beautiDocFor(nexus);
        File out = File.createTempFile("readCountBeautiCompare", ".xml");
        out.deleteOnExit();
        doc.save(out);

        MCMC mcmc = (MCMC) new XMLParser().parseFile(out);
        mcmc.chainLengthInput.setValue(0L, mcmc);
        mcmc.initAndValidate();
        mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());

        ReadCountTreeLikelihood generated = null;
        for (Distribution d : ((CompoundDistribution) mcmc.posteriorInput.get()).pDistributions.get()) {
            generated = findLikelihood(d);
            if (generated != null) {
                break;
            }
        }
        assertTrue(generated != null, "no ReadCountTreeLikelihood in the generated analysis");

        // rebuild the same likelihood from scratch: independently imported data, the template's
        // parameter values written out here rather than read back from the generated model, and
        // the tree BEAUti initialised (the starting tree is random, so it has to be shared)
        ReadCountAlignment data = (ReadCountAlignment) new ReadCountNexusImporter()
                .loadFile(nexus).get(0);

        SimplexParam freqs = new SimplexParam(uniform(16, 0.0625));
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        RealVectorParam nucRates = new RealVectorParam(uniform(6, 1.0 / 6.0), PositiveReal.INSTANCE);
        nucRates.setInputValue("keys", "AC AG AT CG CT GT");
        nucRates.initAndValidate();
        GT16 gt16 = new GT16();
        gt16.initByName("nucRates", nucRates, "frequencies", frequencies);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", gt16, "gammaCategoryCount", 0,
                "mutationRate", scalar(1.0), "shape", scalar(1.0),
                "proportionInvariant", new RealScalarParam(0.0, UnitInterval.INSTANCE));

        LikelihoodReadCountModel rcm = new LikelihoodReadCountModel();
        rcm.initByName("readCount", data,
                "epsilon", new RealScalarParam(0.01, UnitInterval.INSTANCE),
                "delta", new RealScalarParam(0.05, UnitInterval.INSTANCE),
                "t", scalar(5.0), "v", scalar(1.0),
                "s", new RealVectorParam(uniform(data.getTaxaNumber(), 1.0), PositiveReal.INSTANCE),
                "w1", scalar(100.0), "w2", scalar(2.0));

        ReadCountTreeLikelihood handBuilt = new ReadCountTreeLikelihood();
        handBuilt.initByName("data", data, "tree", generated.treeInput.get(),
                "siteModel", siteModel, "readCountModel", rcm);

        assertEquals(generated.calculateLogP(), handBuilt.calculateLogP(), 1e-9,
                "the BEAUti-generated likelihood differs from an equivalent hand-built one");
    }

    private static ReadCountTreeLikelihood findLikelihood(Distribution d) {
        if (d instanceof ReadCountTreeLikelihood) {
            return (ReadCountTreeLikelihood) d;
        }
        if (d instanceof CompoundDistribution) {
            for (Distribution child : ((CompoundDistribution) d).pDistributions.get()) {
                ReadCountTreeLikelihood found = findLikelihood(child);
                if (found != null) {
                    return found;
                }
            }
        }
        return null;
    }

    private static RealScalarParam scalar(double v) {
        return new RealScalarParam(v, PositiveReal.INSTANCE);
    }

    private static double[] uniform(int n, double v) {
        double[] a = new double[n];
        java.util.Arrays.fill(a, v);
        return a;
    }
}
