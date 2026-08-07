package phylonco.beast.evolution.likelihood;


import beagle.BeagleFactory;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.substitutionmodel.Frequencies;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;
import phylonco.beast.evolution.substitutionmodel.GT10;

import java.util.Arrays;

import static beast.pkgmgmt.BEASTClassLoader.addServices;
import static org.junit.jupiter.api.Assertions.*;



/**
 * Validates that {@link ReadCountTreeLikelihood} (which integrates the latent genotypes out via
 * read-count tip partials) equals a hand-computed Felsenstein integral over genotypes that uses
 * the same per-genotype read-count likelihoods. This cross-checks the partials wiring into the
 * Beer likelihood core: genotype state ordering, branch direction, root frequencies as the
 * genotype prior, and the per-(tip,site) offset bookkeeping.
 */
public class BeagleReadCountTreeLikelihoodTest {

    private static boolean beagleFound = false;
    // executed only once, before the first test
    @BeforeAll
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

    private static final double DELTA = 1e-9;
    private static final int N = 10; // GT10 genotype states

    // taxa "a", "b"; tree (a:0.3, b:0.5)
    private static final double BRANCH_A = 0.3;
    private static final double BRANCH_B = 0.5;

    // read-count model parameters
    private static final String EPSILON = "0.06";
    private static final String DELTA_DROPOUT = "0.5";
    private static final String T_PARAM = "9.996182050184155";
    private static final String V_PARAM = "1.0670434040009762";
    private static final double[] S_PARAM = {1.0399635911708527, 1.0419228814287969};
    private static final String W1 = "100.0";
    private static final String W2 = "2.0";

    @BeforeAll
    public static void setUpClass() {
        addServices("version.xml");
    }

    /** Builds the GT10 substitution + site model with uniform frequencies and fixed nuc rates. */
    private SiteModel buildSiteModel() {
        double[] piVals = new double[N];
        Arrays.fill(piVals, 1.0 / N);
        SimplexParam pi = new SimplexParam(piVals);
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", pi, "estimate", false);
        freqs.initAndValidate();

        RealVectorParam nucRates = new RealVectorParam(
                new double[]{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, PositiveReal.INSTANCE);
        nucRates.setInputValue("keys", "AC AG AT CG CT GT");
        nucRates.initAndValidate();

        GT10 subsModel = new GT10();
        subsModel.initByName("nucRates", nucRates, "frequencies", freqs);
        subsModel.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "mutationRate", new RealScalarParam(1.0, PositiveReal.INSTANCE),
                "gammaCategoryCount", 1, "substModel", subsModel);
        siteModel.initAndValidate();
        return siteModel;
    }

    /** Constant scaffold alignment: GT10 datatype, taxa a/b, distinct columns so patterns==sites. */
    private Alignment buildScaffold(int nSites) {
        // "AC" -> sites (AA, CC); "GT" -> sites (GG, TT): distinct columns per site, no compression
        Alignment data = new Alignment();
        data.initByName(
                "sequence", new Sequence("a", "AC".substring(0, nSites)),
                "sequence", new Sequence("b", "GT".substring(0, nSites)),
                "dataType", "nucleotideDiploid10");
        return data;
    }

    /** A positive-real scalar parameter from a decimal string (spec params no longer auto-convert Strings). */
    private static RealScalarParam scalar(String value) {
        return new RealScalarParam(Double.parseDouble(value), PositiveReal.INSTANCE);
    }

    private LikelihoodReadCountModel buildReadCountModel(Alignment scaffold, String readCountStr) {
        ReadCount readCount = new ReadCount();
        readCount.initByName("value", readCountStr, "taxaNames", "a b");

        LikelihoodReadCountModel rcm = new LikelihoodReadCountModel();
        rcm.initByName(
                "alignment", scaffold,
                "readCount", readCount,
                "epsilon", scalar(EPSILON),
                "delta", scalar(DELTA_DROPOUT),
                "t", scalar(T_PARAM),
                "v", scalar(V_PARAM),
                "s", new RealVectorParam(S_PARAM, PositiveReal.INSTANCE),
                "w1", scalar(W1),
                "w2", scalar(W2));
        return rcm;
    }

    private BeagleReadCountTreeLikelihood buildLikelihood(Alignment scaffold, LikelihoodReadCountModel rcm,
                                                    SiteModel siteModel, TreeParser tree) {
        BeagleReadCountTreeLikelihood likelihood = new BeagleReadCountTreeLikelihood();
        likelihood.initByName(
                "data", scaffold,
                "tree", tree,
                "siteModel", siteModel,
                "readCountModel", rcm);
        return likelihood;
    }

    private TreeParser buildTree(Alignment scaffold) {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", scaffold, "newick", "(a:" + BRANCH_A + ", b:" + BRANCH_B + ");",
                "IsLabelledNewick", true);
        return tree;
    }

    /**
     * Hand Felsenstein reference, done per-site in log space with the same max-offset trick the
     * implementation uses, so it is both the ground truth and numerically robust under underflow.
     */
    private double referenceLogP(Alignment scaffold, LikelihoodReadCountModel rcm,
                                 SiteModel siteModel, TreeParser tree, int nSites) {
        rcm.initialize();
        GT10 subs = (GT10) siteModel.substModelInput.get();
        double[] freq = subs.getFrequencies();
        int[] a2rc = rcm.getAlignToRCIndex();

        Node nodeA = nodeOf(tree, "a");
        Node nodeB = nodeOf(tree, "b");
        // use the branch lengths BEAST actually assigned (robust to non-ultrametric height resolution)
        double[] Pa = transitionMatrix(subs, nodeA, nodeA.getLength());
        double[] Pb = transitionMatrix(subs, nodeB, nodeB.getLength());

        int alignA = scaffold.getTaxonIndex("a");
        int alignB = scaffold.getTaxonIndex("b");

        double logRef = 0.0;
        for (int j = 0; j < nSites; j++) {
            double[] logA = perGenotypeLog(rcm, alignA, a2rc[alignA], j);
            double[] logB = perGenotypeLog(rcm, alignB, a2rc[alignB], j);
            double offA = max(logA);
            double offB = max(logB);
            double[] La = scaledExp(logA, offA);
            double[] Lb = scaledExp(logB, offB);

            double Lj = 0.0;
            for (int r = 0; r < N; r++) {
                double sa = 0.0, sb = 0.0;
                for (int g = 0; g < N; g++) {
                    sa += Pa[r * N + g] * La[g];
                    sb += Pb[r * N + g] * Lb[g];
                }
                Lj += freq[r] * sa * sb;
            }
            logRef += Math.log(Lj) + offA + offB;
        }
        return logRef;
    }

    private double[] perGenotypeLog(LikelihoodReadCountModel rcm, int alignIdx, int rcIdx, int site) {
        int[] reads = rcm.getReadCount().getReadCounts(rcIdx, site);
        int coverage = rcm.getCoverage(alignIdx, site);
        double[] out = new double[N];
        for (int g = 0; g < N; g++) out[g] = rcm.logLiklihoodRC(g, reads, coverage, rcIdx);
        return out;
    }

    private double[] transitionMatrix(GT10 subs, Node node, double branchLen) {
        double[] mat = new double[N * N];
        subs.getTransitionProbabilities(node, branchLen, 0.0, 1.0, mat);
        return mat;
    }

    private Node nodeOf(TreeParser tree, String id) {
        for (Node n : tree.getExternalNodes()) if (id.equals(n.getID())) return n;
        throw new RuntimeException("no node " + id);
    }

    private static double max(double[] x) {
        double m = Double.NEGATIVE_INFINITY;
        for (double v : x) if (v > m) m = v;
        return m;
    }

    private static double[] scaledExp(double[] log, double off) {
        double[] out = new double[log.length];
        for (int i = 0; i < log.length; i++) out[i] = Math.exp(log[i] - off);
        return out;
    }

    /** Core: 2 taxa x 2 sites, GT10. Integrated likelihood == hand Felsenstein integral. */
    @Test
    public void testIntegratedEqualsBruteForceTwoSite() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            // TODO: Kylie to check - silently returning here reports this test as "passed" without running any
            // assertions when BEAGLE isn't installed, masking a missing/broken BEAGLE setup as
            // green. Check whether this should fail (or use Assumptions.assumeTrue, which reports
            // as "skipped" rather than "passed") instead.
            return;
        }
        int nSites = 2;
        // a: site0 / site1 ; b: site0 / site1   (A:C:G:T per site, sites separated by ",")
        String readCountStr = "3:0:1:8,0:12:2:0\n7:1:0:4,0:1:9:1";
        Alignment scaffold = buildScaffold(nSites);
        LikelihoodReadCountModel rcm = buildReadCountModel(scaffold, readCountStr);
        SiteModel siteModel = buildSiteModel();
        TreeParser tree = buildTree(scaffold);

        BeagleReadCountTreeLikelihood likelihood = buildLikelihood(scaffold, rcm, siteModel, tree);

        double logP = likelihood.calculateLogP();
        double logRef = referenceLogP(scaffold, rcm, siteModel, tree, nSites);
        assertEquals(logRef, logP, DELTA);
    }

    /** Underflow: a very-high-coverage site makes raw read-count likelihoods strongly negative;
     *  the integrated result must still match the (offset-robust) reference exactly. */
    @Test
    public void testUnderflowHighCoverage() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            // TODO: Kylie to check - silently returning here reports this test as "passed" without running any
            // assertions when BEAGLE isn't installed, masking a missing/broken BEAGLE setup as
            // green. Check whether this should fail (or use Assumptions.assumeTrue, which reports
            // as "skipped" rather than "passed") instead.
            return;
        }
        int nSites = 1;
        // coverage ~300 at one site for each cell
        String readCountStr = "120:5:170:5\n3:150:7:140";
        Alignment scaffold = buildScaffold(nSites);
        LikelihoodReadCountModel rcm = buildReadCountModel(scaffold, readCountStr);
        SiteModel siteModel = buildSiteModel();
        TreeParser tree = buildTree(scaffold);

        BeagleReadCountTreeLikelihood likelihood = buildLikelihood(scaffold, rcm, siteModel, tree);

        double logP = likelihood.calculateLogP();
        double logRef = referenceLogP(scaffold, rcm, siteModel, tree, nSites);
        assertEquals(logRef, logP, DELTA);
        // sanity: a real (finite) number, not -Inf/NaN
        assertTrue(Double.isFinite(logP));
    }

    /** Stale-cache guard (must-fix #1): perturbing a read-count param then restoring it must
     *  return exactly the original logP — proves initialize() is refreshed on every rebuild. */
    @Test
    public void testParamRoundTripConsistency() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            // TODO: Kylie to check - silently returning here reports this test as "passed" without running any
            // assertions when BEAGLE isn't installed, masking a missing/broken BEAGLE setup as
            // green. Check whether this should fail (or use Assumptions.assumeTrue, which reports
            // as "skipped" rather than "passed") instead.
            return;
        }
        int nSites = 2;
        String readCountStr = "3:0:1:8,0:12:2:0\n7:1:0:4,0:1:9:1";
        Alignment scaffold = buildScaffold(nSites);
        LikelihoodReadCountModel rcm = buildReadCountModel(scaffold, readCountStr);
        SiteModel siteModel = buildSiteModel();
        TreeParser tree = buildTree(scaffold);

        BeagleReadCountTreeLikelihood likelihood = buildLikelihood(scaffold, rcm, siteModel, tree);
        RealScalarParam epsilon = (RealScalarParam) rcm.epsilonInput.get();

        double logP0 = likelihood.calculateLogP();

        epsilon.set(0.2);
        rcm.requiresRecalculation();
        likelihood.requiresRecalculation();
        double logP1 = likelihood.calculateLogP();
        //System.out.println("1. Expected: " + logP0 + ", Actual: "+ logP1);
        System.out.println("logP0: " + logP0 );
        System.out.println("logP1: " + logP0 );
        assertNotEquals(logP0, logP1, DELTA); // the move actually changed something

        epsilon.set(Double.parseDouble(EPSILON));
        rcm.requiresRecalculation();
        likelihood.requiresRecalculation();
        double logP2 = likelihood.calculateLogP();
        //System.out.println("2. Expected: " + logP0 + ", Actual: "+ logP2);
        assertEquals(logP0, logP2, DELTA);
    }

    /** store()/restore() of globalLogOffset (must-fix #3): after store, a param change, then
     *  restore (with the param reverted), a fresh logP must equal the pre-move value. */
    @Test
    public void testStoreRestoreOffset() {
        if (beagleFound == false) {
            System.err.println("Warning: Cannot find beagle resources! Skipping test.");
            // TODO: Kylie to check - silently returning here reports this test as "passed" without running any
            // assertions when BEAGLE isn't installed, masking a missing/broken BEAGLE setup as
            // green. Check whether this should fail (or use Assumptions.assumeTrue, which reports
            // as "skipped" rather than "passed") instead.
            return;
        }
        int nSites = 2;
        String readCountStr = "3:0:1:8,0:12:2:0\n7:1:0:4,0:1:9:1";
        Alignment scaffold = buildScaffold(nSites);
        LikelihoodReadCountModel rcm = buildReadCountModel(scaffold, readCountStr);
        SiteModel siteModel = buildSiteModel();
        TreeParser tree = buildTree(scaffold);

        BeagleReadCountTreeLikelihood likelihood = buildLikelihood(scaffold, rcm, siteModel, tree);
        RealScalarParam epsilon = (RealScalarParam) rcm.epsilonInput.get();

        double logP0 = likelihood.calculateLogP();

        likelihood.store();
        rcm.store();
        epsilon.set(0.2);
        likelihood.requiresRecalculation();
        likelihood.calculateLogP();

        // reject: revert param and restore the calculation nodes
        epsilon.set(Double.parseDouble(EPSILON));
        likelihood.restore();
        rcm.restore();
        double logP2 = likelihood.calculateLogP();
        assertEquals(logP0, logP2, DELTA);
    }
}
