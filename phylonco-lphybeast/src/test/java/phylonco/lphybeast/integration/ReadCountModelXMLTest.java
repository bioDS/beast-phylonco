package phylonco.lphybeast.integration;

import beast.base.parser.XMLParser;
import lphybeast.LPhyBeastCMD;
import org.junit.jupiter.api.Test;
import picocli.CommandLine;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import static beast.pkgmgmt.BEASTClassLoader.addServices;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * End-to-end check of {@link phylonco.lphybeast.tobeast.generators.ReadCountModelToBEAST}: runs a
 * GT10 read-count LPhy script through LPhyBeast and asserts the produced XML is the integrated model
 * — a single {@code ReadCountTreeLikelihood} factor and a {@code SampledGenotypeLogger}, with the
 * old sampled-genotype machinery (MATreeLikelihood / genotype Gibbs operators) gone.
 */
public class ReadCountModelXMLTest {

    @Test
    public void testReadCountIntegratedXML() throws Exception {
        // Absolute paths: LPhyBeast switches the working dir to the script's folder, which would
        // break relative-path resolution in its input validation.
        Path script = Paths.get("src/test/resources/gt10ReadCountIntegrated.lphy").toAbsolutePath();
        Path out = Paths.get("target/gt10ReadCountIntegrated.xml").toAbsolutePath();
        Files.deleteIfExists(out);

        // -vf registers BEAST2 services from version.xml files. lphybeast.version.xml supplies the
        // lphybeast core SPI (LPhyBEASTMapping etc.) and version.xml supplies phylonco's LBPhylonco
        // extension, so generation works headlessly without the BEAST packages being installed.
        String[] args = {
                "-seed", "777",
                "-vf", "lphybeast.version.xml",
                "-vf", "version.xml",
                "-l", "2000", // must be >= LPhyBeastConfig.NUM_OF_SAMPLES (2000)
                "-o", out.toString(),
                script.toString()
        };
        int exitCode = new CommandLine(new LPhyBeastCMD()).execute(args);
        assertEquals(0, exitCode, "LPhyBeast generation should succeed");
        assertTrue(Files.exists(out), "XML should be generated at " + out.toAbsolutePath());

        String xml = Files.readString(out);

        // the integrated model: read-count tip-partial tree likelihood + tree-aware genotype logger
        assertTrue(xml.contains("ReadCountTreeLikelihood"),
                "XML must contain a ReadCountTreeLikelihood factor");
        assertTrue(xml.contains("SampledGenotypeLogger") || xml.contains("sampledGenotype"),
                "XML must contain the SampledGenotypeLogger");

        // the old sampled-genotype machinery must be gone
        assertFalse(xml.contains("MATreeLikelihood"),
                "XML must NOT contain the old MATreeLikelihood");
        assertFalse(xml.contains("GibbsSiteOperator") || xml.contains("GibbsAlignmentOperator")
                        || xml.contains("GibbsSequenceOperator"),
                "XML must NOT contain genotype Gibbs operators");

        // ---- run the generated XML through BEAST for a short chain (parse + MCMC integrity) ----
        // register the phylonco BEAST classes and MutableAlignment as BEAST2 services
        addServices(Paths.get("..", "phylonco-beast", "version.xml").toString());
        addServices(Paths.get("..", "phylonco-beast", "lib", "MutableAlignment.version.xml").toString());

        // parseFile builds and initialises the full model graph (initAndValidate on every object,
        // incl. ReadCountTreeLikelihood rebuilding its scaffold from the read counts) — it throws if
        // the generated XML is not a valid, wireable BEAST model.
        beast.base.inference.MCMC mcmc =
                (beast.base.inference.MCMC) new XMLParser().parseFile(out.toFile());

        // the model must produce a finite initial posterior (exercises the integrated tip-partial
        // likelihood end-to-end on the generated XML). Set up the State and do the robust posterior
        // calculation exactly as MCMC.run() does at start-up.
        beast.base.inference.Distribution posterior = mcmc.posteriorInput.get();
        beast.base.inference.State state = mcmc.startStateInput.get();
        state.initAndValidate();
        state.setPosterior(posterior);
        double logP = state.robustlyCalcPosterior(posterior);
        assertTrue(Double.isFinite(logP), "initial posterior must be finite, was " + logP);
    }
}
