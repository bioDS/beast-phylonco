package phylonco.beast.evolution.likelihood;


import beagle.BeagleFactory;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.Binary;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.substitutionmodel.Frequencies;
import beast.base.spec.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import org.junit.BeforeClass;
import org.junit.Test;
import phylonco.beast.evolution.datatype.NucleotideDiploid16;
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
    private static boolean beagleFound = false;

    // executed only once, before the first test
    @BeforeClass
    public static void setUpBeagle() {
        // add -Djava.library.path="$LD_LIBRARY_PATH:/usr/local/lib" before running tests
        // "/usr/local/lib" is the Beagle lib location in this case
        System.out.println("Warning: Add -Djava.library.path=YOUR_LIB_PATH before running tests!");
        System.out.println("java.library.path = " + System.getProperty("java.library.path"));
        System.out.println(BeagleFactory.getVersionInformation());
        if (BeagleFactory.getResourceDetails().isEmpty()) {
            beagleFound = false;
            System.err.println("Warning: Cannot find beagle resources! Skipping BeagleTreeLikelihoodWithErrorTest");
        } else {
            beagleFound = true;
            // print Beagle information
            System.out.println("\n--- BEAGLE RESOURCES ---\n");
            for (beagle.ResourceDetails details : BeagleFactory.getResourceDetails()) {
                System.out.println(details.toString());
            }
        }
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
        siteModel.initByName(
                "mutationRate", new RealScalarParam(1.0, PositiveReal.INSTANCE),
                "gammaCategoryCount", 1,
                "substModel", subsModel);
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
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            return;
        }
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
        errorModel.initByName(
                "epsilon", new RealScalarParam(0.1, UnitInterval.INSTANCE),
                "datatype", datatype);
        errorModel.initAndValidate();

        double logP = getLogLikelihood(data, tree, siteModel, errorModel);
        double expectedLogP = -2.3063595712034233;
        assertEquals(expectedLogP, logP, DELTA);
    }

    private double calculateLikelihoodBinary(String seq, double alpha, double beta) {
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
        subsModel.initByName("lambda", new RealScalarParam(2.0, PositiveReal.INSTANCE));
        subsModel.initAndValidate();

        SiteModel siteModel = getSiteModel(subsModel);

        Binary datatype = new Binary();

        BinaryErrorModel errorModel = new BinaryErrorModel();
        errorModel.initByName(
                "alpha", new RealScalarParam(alpha, UnitInterval.INSTANCE),
                "beta", new RealScalarParam(beta, UnitInterval.INSTANCE),
                "datatype", datatype);
        errorModel.initAndValidate();

        return getLogLikelihood(data, tree, siteModel, errorModel);
    }

    @Test
    public void testBinaryLikelihoodSmallNoError() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            return;
        }
        double logP = calculateLikelihoodBinary("00", 0.0, 0.0);
        double expectedLogP = -0.7595722922504291;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testBinaryLikelihoodSmallWithErrorCase0() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            return;
        }
        double logP = calculateLikelihoodBinary("00", 0.1, 0.2);
        double expectedLogP = -0.78543518416993563;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testBinaryLikelihoodSmallWithErrorCase1() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            return;
        }
        double logP = calculateLikelihoodBinary("11", 0.1, 0.2);
        double expectedLogP = -2.0989268283365146;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testBinaryLikelihoodSmallWithErrorCase2() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            return;
        }
        double logP = calculateLikelihoodBinary("01", 0.1, 0.2);
        double expectedLogP = -1.5571044248279775;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testBinaryLikelihoodSmallTotalProbability() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            return;
        }
        double logP1 = calculateLikelihoodBinary("00", 0.1, 0.2);
        double logP2 = calculateLikelihoodBinary("01", 0.1, 0.2);
        double logP3 = calculateLikelihoodBinary("10", 0.1, 0.2);
        double logP4 = calculateLikelihoodBinary("11", 0.1, 0.2);
        double probSum = Math.exp(logP1) + Math.exp(logP2) + Math.exp(logP3) + Math.exp(logP4);
        assertEquals(1.0, probSum, DELTA);
    }

    private double calculateLikelihoodGT16(String seq, double epsilon, double delta) {
        Alignment data = new Alignment();
        Sequence seqA = new Sequence("a", seq.substring(0, 1));
        Sequence seqB = new Sequence("b", seq.substring(1));
        data.initByName(
                "sequence", seqA,
                "sequence", seqB,
                "dataType", "nucleotideDiploid16"
        );

        TreeParser tree = getTree(data);

        double[] pi = new double[16];
        Arrays.fill(pi, 1.0 / 16);
        SimplexParam f = new SimplexParam(pi);
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f, "estimate", false);
        freqs.initAndValidate();

        double[] rates = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        RealVectorParam nucRates = new RealVectorParam(rates, PositiveReal.INSTANCE);
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
        errorModel.initByName(
                "epsilon", new RealScalarParam(epsilon, UnitInterval.INSTANCE),
                "delta", new RealScalarParam(delta, UnitInterval.INSTANCE),
                "datatype", datatype);
        errorModel.initAndValidate();

        return getLogLikelihood(data, tree, siteModel, errorModel);
    }

    @Test
    public void testGT16ErrorLikelihoodCase0() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            return;
        }
        double logP = calculateLikelihoodGT16("00", 0.1, 0.2);
        double expectedLogP = -3.2683402019565975;
        assertEquals(expectedLogP, logP, DELTA);
    }

    @Test
    public void testGT16ErrorLikelihoodCase1() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            return;
        }
        double logP = calculateLikelihoodGT16("01", 0.1, 0.2);
        double expectedLogP = -5.1071258693509041;
        assertEquals(expectedLogP, logP, DELTA);
    }
}