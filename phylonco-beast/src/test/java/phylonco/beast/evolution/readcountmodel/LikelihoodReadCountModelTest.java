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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.junit.Assert.assertEquals;


public class LikelihoodReadCountModelTest {

    final double DELTA = 1e-5;

    @BeforeEach
    public void setUp() {
        String versionFile = "version.xml";
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
        LikelihoodReadCountModel likelihoodReadCountModel = new LikelihoodReadCountModel();

        // read from file
        Double epsilon = 0.06;
        Double delta = 0.5;
        Double t = 9.996182050184155;
        Double v = 1.0670434040009762;
        Double[] s = new Double[]{1.0399635911708527, 1.0419228814287969};
        Double w1 = 10.0;
        Double w2 = 2.0;

        Path dir = Path.of("src","test", "resources");
        File alignmentFile = Paths.get(dir.toString(),"gt16ReadCountModel_A.nexus").toFile();
        File readCountFile = Paths.get(dir.toString(),"readCountNumbers.txt").toFile();

        Alignment alignment = getAlignment(alignmentFile);
        ReadCount readCounts = getReadCounts(readCountFile);

        // real parameter array
        RealParameter sParam = new RealParameter(s);

        // init params
        likelihoodReadCountModel.setInputValue("alignment", alignment);
        likelihoodReadCountModel.setInputValue("readCount", readCounts);
        likelihoodReadCountModel.setInputValue("epsilon", epsilon.toString());
        likelihoodReadCountModel.setInputValue("delta", delta.toString());
        likelihoodReadCountModel.setInputValue("t", t.toString());
        likelihoodReadCountModel.setInputValue("v", v.toString());
        likelihoodReadCountModel.setInputValue("s", sParam);
        likelihoodReadCountModel.setInputValue("w1", w1.toString());
        likelihoodReadCountModel.setInputValue("w2", w2.toString());

        // ...

        likelihoodReadCountModel.initAndValidate();

        double observedLogP = likelihoodReadCountModel.calculateLogP();
        //System.out.println(observedLogP);
        double expectedLogP = -47.1349957401387;

        assertEquals(expectedLogP, observedLogP, DELTA);

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

    private Alignment getAlignment(File file) {
        System.out.println("Processing " + file);
        NexusParser parser = new NexusParser();
        try {
            parser.parseFile(file);
            System.out.println("Done " + file);
            return parser.m_alignment;
        } catch (Exception e) {
            throw new RuntimeException("ExampleNexusParsing::Failed for " + file
                    + ": " + e.getMessage());
        }
    }

    private ReadCount getReadCounts(File file) throws IOException {
        ReadCount readCount;
        int ntaxa;
        int nchar;
        int lineCount = 0;
        ArrayList<Integer> numbers = new ArrayList<>();
        // File reader
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        while ((line = reader.readLine()) != null) {
            Pattern pattern = Pattern.compile("\\d+");
            Matcher matcher = pattern.matcher(line);
            while (matcher.find()) {
                numbers.add(Integer.parseInt(matcher.group()));
            }
            lineCount++;
        }
        int[] numArray = numbers.stream().mapToInt(i -> i).toArray();
        ntaxa = lineCount;
        nchar = numbers.size()/lineCount/4;
        readCount = new ReadCount(ntaxa, nchar);
        int[] Count = new int[4];
        for (int i = 0; i < ntaxa; i++) {
            for (int j = 0; j < nchar; j++) {
                for (int k = 0; k < 4; k++) {
                    Count[k] = numArray[i*nchar*4+j*4+k];
                }
                readCount.setReadCounts(i, j, Count);
            }
        }

        return readCount;
    }

}
