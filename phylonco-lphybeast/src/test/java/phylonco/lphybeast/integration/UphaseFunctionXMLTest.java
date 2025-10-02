package phylonco.lphybeast.integration;

import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.function.io.ReadNexus;
import lphy.core.model.Value;
import lphybeast.LPhyBeastCMD;
import org.junit.jupiter.api.Test;
import picocli.CommandLine;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class UphaseFunctionXMLTest {


    @Test
    public void test_UnphaseXMLExample() {
        try {
            String[] args = {
                    "-seed", "666",
                    "-vf", "version.xml",
                    "-r", "1",
                    "src/test/resources/gt16CoalErrModelUnphase.lphy"
                    };
            int exitCode = new CommandLine(new LPhyBeastCMD()).execute(args);
            Value<String> expectedN = new Value<>("expectedNexus","exampleGt16CoalErrModelUnphase_true_D.nexus");
            Value<String> generatedN = new Value<>("generatedNexus","gt16CoalErrModelUnphase_true_D.nexus");
            String expected = new ReadNexus(expectedN,null).apply().value().toJSON();
            String generated = new ReadNexus(generatedN,null).apply().value().toJSON();
            assertEquals(true, expected.equals(generated));
        } catch (Exception e) {
            System.out.println("exception thrown ");
            System.out.println(e.getMessage());
        }
    }




}
