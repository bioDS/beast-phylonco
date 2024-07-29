package phylonco.beast.evolution.readcountmodel;

import beast.base.evolution.alignment.Alignment;
import beast.base.inference.parameter.RealParameter;
import beast.base.parser.NexusParser;
import beast.pkgmgmt.BEASTClassLoader;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import phylonco.beast.evolution.datatype.ReadCount;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class ReadCountModelTest {

    final double DELTA = 1e-6;

    @BeforeEach
    public void setUp() {
        String versionFile = "../phylonco-lphybeast/version.xml";
        BEASTClassLoader.addServices(versionFile);
    }

    @Test
    public void testAlignment() {
        Alignment alignment = new Alignment();
    }


    /**
     * Test log likelihood of read count model using data generated usign LPhy script:
     *
     *
     * Expected log likelihoood calculated using R script: calcLogP.R
     */
    @Test
    public void testReadCountModel() throws IOException {
        ReadCountModel readCountModel = new ReadCountModel();

        // read from file
        Double epsilon = 0.1;
        Double delta = 0.24;
        Double t = 10.533354280335304;
        Double v = 1.0158887508905496;
        Double[] s = new Double[]{1.040558132021848, 1.0411303267252416, 1.0404430631662909, 1.0403929520996686};
        Double w = 100.0;


        String alignmentFile = "/Users/yxia415/Desktop/data/alignment.nexus";
        String readCountFile = "/Users/yxia415/Desktop/data/readCount.txt";
        Alignment alignment = getAlignment(alignmentFile);
        ReadCount readCounts = getReadCounts(readCountFile);

        // real parameter array
        RealParameter sParam = new RealParameter(s);

        // init params
        readCountModel.setInputValue("alignment", alignment);
        readCountModel.setInputValue("readCount", readCounts);
        readCountModel.setInputValue("epsilon", epsilon.toString());
        readCountModel.setInputValue("delta", delta.toString());
        readCountModel.setInputValue("t", t.toString());
        readCountModel.setInputValue("v", v.toString());
        readCountModel.setInputValue("s", sParam);
        readCountModel.setInputValue("w", w.toString());

        // ...

        readCountModel.initAndValidate();

        double observedLogP = readCountModel.calculateLogP();
        double expectedLogP = -1.0;

//        assertEquals(expectedLogP, observedLogP, DELTA);

//        public Input<Alignment> alignmentInput = new Input<>("alignment", "alignment");
//        public Input<ReadCount> readCountInput = new Input<>("readCount", "nucleotide read counts");
//
//        // epsilon, allelic dropout, ... parameters
//        public Input<RealParameter> epsilonInput = new Input<>("epsilon", "sequencing error");
//        public Input<RealParameter> deltaInput = new Input<>("delta", "allelic dropout probability");
//        public Input<RealParameter> tInput = new Input<>("t", "mean of allelic coverage");
//        public Input<RealParameter> vInput = new Input<>("v", "variance of allelic coverage");
//        public Input<RealParameter> sInput = new Input<>("s", "size factor of cell");
//        public Input<RealParameter> wInput = new Input<>("w", "overdispersion parameter of Dirichlet multinomial distribution");


    }

    private Alignment getAlignment(String fileName) {
        System.out.println("Processing " + fileName);
        NexusParser parser = new NexusParser();
        try {
            parser.parseFile(new File(fileName));
            return parser.m_alignment;
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("ExampleNexusParsing::Failed for " + fileName
                    + ": " + e.getMessage());
        }
        System.out.println("Done " + fileName);
        return null;
    }

    private ReadCount getReadCounts(String fileName) throws IOException {
        // File reader
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        String line = reader.readLine();
        System.out.println(line);
        // for loop for reading each line
        // split by tabs for each site
        // split by commas for each nucleotide
        // put read count numbers into int array a
        int[] a; // hard code if you want??
        // ReadCount r = new ReadCount(a);
        return null;
    }

}
