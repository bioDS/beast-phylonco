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

public class ExchangeGibbsOperatorShortTest {

    @BeforeEach
    public void setUp() {
        addServices("version.xml");
        addServices(Paths.get("lib", "MutableAlignment.version.xml").toString());
    }

    @Test
    public void test() throws Exception {
        File xmlFile = Paths.get("src","test","resources","ExchangeGibbsOperatorTest_short.xml").toFile();
        BeastMCMC beastMCMC = new BeastMCMC();
        List<String> args = new ArrayList<>();
        args.add("-overwrite");
        args.add(xmlFile.getAbsolutePath());
        beastMCMC.parseArgs(args.toArray(new String[0]));
        beastMCMC.run();
    }
}
