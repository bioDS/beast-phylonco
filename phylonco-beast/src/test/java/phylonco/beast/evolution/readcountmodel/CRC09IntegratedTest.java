package phylonco.beast.evolution.readcountmodel;

import beastfx.app.beast.BeastMCMC;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static beast.pkgmgmt.BEASTClassLoader.addServices;

/**
 * End-to-end check of the integrated-genotype model on real CRC09 data: genotypes are integrated
 * out via {@link phylonco.beast.evolution.likelihood.ReadCountTreeLikelihood} tip partials, with no
 * MutableAlignment StateNode and no Gibbs operators. A short chain validates that the model parses,
 * initialises (state-count / identity-pattern / proportionInvariant asserts pass), and runs the MCMC
 * without tripping the BEAST likelihood-integrity check.
 */
public class CRC09IntegratedTest {

    @BeforeEach
    public void setUp() {
        addServices("version.xml");
        addServices(Paths.get("lib", "MutableAlignment.version.xml").toString());
    }

    @Test
    public void test() throws Exception {
        File xmlFile = Paths.get("src", "test", "resources", "CRC09_integrated.xml").toFile();
        BeastMCMC beastMCMC = new BeastMCMC();
        List<String> args = new ArrayList<>();
        args.add("-overwrite");
        args.add(xmlFile.getAbsolutePath());
        beastMCMC.parseArgs(args.toArray(new String[0]));
        beastMCMC.run();
    }
}
