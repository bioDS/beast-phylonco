package beast.evolution.likelihood;


import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.Binary;
import beast.base.evolution.datatype.Nucleotide;
import beast.evolution.datatype.NucleotideDiploid16;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.TreeParser;
import org.junit.BeforeClass;
import org.junit.Test;
import phylonco.beast.evolution.errormodel.BinaryErrorModel;
import phylonco.beast.evolution.errormodel.ErrorModel;
import phylonco.beast.evolution.errormodel.ErrorModelBase;
import phylonco.beast.evolution.errormodel.GT16ErrorModel;
import phylonco.beast.evolution.substitutionmodel.BinarySubstitutionModel;

import java.util.Arrays;

import static org.junit.Assert.*;

/**
 * This requires Beagle installed, and the Beagle lib location,
 * before running the tests.
 * E.g. <code>-Djava.library.path="$LD_LIBRARY_PATH:/usr/local/lib"</code>
 * @author Walter Xie
 */
public class BeagleTreeLikelihoodWithErrorTest {

    private static final double DELTA = 1e-10;
    private static final boolean useGPU = true;

    //executed only once, before the first test
    @BeforeClass
    public static void setUpBeagle() {
        // add -Djava.library.path="$LD_LIBRARY_PATH:/usr/local/lib" before running tests
        // "/usr/local/lib" is the Beagle lib location in this case
        System.out.println("java.library.path = " + System.getProperty("java.library.path"));

        System.out.println(BeagleFactory.getVersionInformation());
        // check if beagle resources cannot be found here
        assertTrue("Cannot find beagle resources !",
                BeagleFactory.getResourceDetails().size() > 0);

        System.out.println("\n--- BEAGLE RESOURCES ---\n");
        for (beagle.ResourceDetails details : BeagleFactory.getResourceDetails())
            System.out.println(details.toString());

        // CPU SSE
        long beagleFlags = BeagleFlag.PROCESSOR_GPU.getMask() | BeagleFlag.VECTOR_SSE.getMask();
        // GPU
        if (useGPU) beagleFlags = BeagleFlag.PROCESSOR_GPU.getMask();
        System.setProperty("beagle.preferred.flags", Long.toString(beagleFlags));
    }

    private static TreeParser getTree(Alignment data) {
        TreeParser tree = new TreeParser();
        tree.initByName(
                "taxa", data,
                "newick", "(a: 0.5, b: 0.5);",
                "IsLabelledNewick", true
        );
        return tree;
    }

    private static SiteModel getSiteModel(SubstitutionModel subsModel) {
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", subsModel);
        siteModel.initAndValidate();
        return siteModel;
    }


    private static double getLogLikelihood(Alignment data, TreeParser tree,
                                           SiteModel siteModel, ErrorModel errorModel) {
        BeagleTreeLikelihoodWithError likelihood = new BeagleTreeLikelihoodWithError();
        likelihood.initByName(
                "data", data,
                "tree", tree,
                "siteModel", siteModel,
                "useAmbiguities", true,
                "useTipLikelihoods", true,
                "errorModel", errorModel);
        assertNotNull("BeagleTreeLikelihoodWithError beagle = " + likelihood.beagle, likelihood.beagle);

        return likelihood.calculateLogP();
    }


    @Test
    public void testJCLikelihoodSmallWithError() {
        Alignment data = new Alignment();
        Sequence seqA = new Sequence("a", "A");
        Sequence seqB = new Sequence("b", "A");
        data.initByName(
                "sequence", seqA,
                "sequence", seqB,
                "dataType", "nucleotide"
        );

        TreeParser tree = getTree(data);

        JukesCantor subsModel = new JukesCantor();
        subsModel.initAndValidate();

        SiteModel siteModel = getSiteModel(subsModel);

        Nucleotide datatype = new Nucleotide();

        ErrorModelBase errorModel = new ErrorModelBase();
        errorModel.initByName("epsilon", "0.1", "datatype", datatype);
        errorModel.initAndValidate();

        double logP = getLogLikelihood(data, tree, siteModel, errorModel);
        double expectedLogP = -2.3063595712034233;
        assertEquals(expectedLogP, logP, DELTA);
    }

    private double calculateLikelihoodBinary(String seq, String alpha, String beta) {
        Alignment data = new Alignment();
        Sequence seqA = new Sequence("a", seq.substring(0, 1));
        Sequence seqB = new Sequence("b", seq.substring(1));
        data.initByName(
                "sequence", seqA,
                "sequence", seqB,
                "dataType", "binary"
        );

        TreeParser tree = getTree(data);

        phylonco.beast.evolution.substitutionmodel.BinarySubstitutionModel subsModel = new BinarySubstitutionModel();
        subsModel.initByName("lambda", "2.0");
        subsModel.initAndValidate();

        SiteModel siteModel = getSiteModel(subsModel);

        Binary datatype = new Binary();

        BinaryErrorModel errorModel = new BinaryErrorModel();
        errorModel.initByName("alpha", alpha, "beta", beta, "datatype", datatype);
        errorModel.initAndValidate();

        return getLogLikelihood(data, tree, siteModel, errorModel);
    }

    @Test
    public void testBinaryLikelihoodSmallNoError() {
        double logP = calculateLikelihoodBinary("00", "0.0", "0.0");
        double expectedLogP = -0.7595722922504291;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testBinaryLikelihoodSmallWithErrorCase0() {
        double logP = calculateLikelihoodBinary("00", "0.1", "0.2");
        double expectedLogP = -0.78543518416993563;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testBinaryLikelihoodSmallWithErrorCase1() {
        double logP = calculateLikelihoodBinary("11", "0.1", "0.2");
        double expectedLogP = -2.0989268283365146;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testBinaryLikelihoodSmallWithErrorCase2() {
        double logP = calculateLikelihoodBinary("01", "0.1", "0.2");
        double expectedLogP = -1.5571044248279775;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testBinaryLikelihoodSmallTotalProbability() {
        double logP1 = calculateLikelihoodBinary("00", "0.1", "0.2");
        double logP2 = calculateLikelihoodBinary("01", "0.1", "0.2");
        double logP3 = calculateLikelihoodBinary("10", "0.1", "0.2");
        double logP4 = calculateLikelihoodBinary("11", "0.1", "0.2");
        double probSum = Math.exp(logP1) + Math.exp(logP2) + Math.exp(logP3) + Math.exp(logP4);
        assertEquals(1.0, probSum, DELTA);
    }

    private double calculateLikelihoodGT16(String seq, String epsilon, String delta) {
        Alignment data = new Alignment();
        Sequence seqA = new Sequence("a", seq.substring(0, 1));
        Sequence seqB = new Sequence("b", seq.substring(1));
        data.initByName(
                "sequence", seqA,
                "sequence", seqB,
                "dataType", "nucleotideDiploid16"
        );

        TreeParser tree = getTree(data);

        Double[] pi = new Double[16];
        Arrays.fill(pi, 1.0 / 16);
        RealParameter f = new RealParameter(pi);
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f, "estimate", false);
        freqs.initAndValidate();

        Double[] rates = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        RealParameter nucRates = new RealParameter(rates);
        nucRates.setInputValue("keys", "AC AG AT CG CT GT");
        nucRates.initAndValidate();

        phylonco.beast.evolution.substitutionmodel.GT16 subsModel = new phylonco.beast.evolution.substitutionmodel.GT16();
        subsModel.initByName(
                "nucRates", nucRates,
                "frequencies", freqs
        );
        subsModel.initAndValidate();

        SiteModel siteModel = getSiteModel(subsModel);

        NucleotideDiploid16 datatype = new NucleotideDiploid16();

        GT16ErrorModel errorModel = new GT16ErrorModel();
        errorModel.initByName("epsilon", epsilon, "delta", delta, "datatype", datatype);
        errorModel.initAndValidate();

        return getLogLikelihood(data, tree, siteModel, errorModel);
    }

    @Test
    public void testGT16ErrorLikelihoodCase0() {
        double logP = calculateLikelihoodGT16("00", "0.1", "0.2");
        double expectedLogP = -3.2683402019565975;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testGT16ErrorLikelihoodCase1() {
        double logP = calculateLikelihoodGT16("01", "0.1", "0.2");
        double expectedLogP = -5.1071258693509041;
        assertEquals(expectedLogP, logP, DELTA);
    }
}