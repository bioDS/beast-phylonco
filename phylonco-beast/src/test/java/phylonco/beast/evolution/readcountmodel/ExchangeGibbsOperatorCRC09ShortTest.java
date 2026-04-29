package phylonco.beast.evolution.readcountmodel;

import beastfx.app.beast.BeastMCMC;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static beast.pkgmgmt.BEASTClassLoader.addServices;

public class ExchangeGibbsOperatorCRC09ShortTest {

    @BeforeEach
    public void setUp() {
        addServices("version.xml");
        addServices(Paths.get("lib", "MutableAlignment.version.xml").toString());
    }

    @Test
    public void test() throws Exception {
        Path dir = Path.of("src", "test", "resources");
        File xmlFile = Paths.get(dir.toString(), "CRC09_short.xml").toFile();
        final BeastMCMC beastMCMC = new BeastMCMC();
        final List<String> MCMCargs = new ArrayList<>();
        MCMCargs.add("-overwrite");
        MCMCargs.add(xmlFile.getAbsolutePath());
        beastMCMC.parseArgs(MCMCargs.toArray(new String[0]));
        beastMCMC.run();
    }
}
