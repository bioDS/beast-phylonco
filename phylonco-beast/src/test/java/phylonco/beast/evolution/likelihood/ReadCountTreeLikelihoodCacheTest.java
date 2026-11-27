package phylonco.beast.evolution.likelihood;

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

import static beast.pkgmgmt.BEASTClassLoader.addServices;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Checks the cached tip partials of {@link ReadCountTreeLikelihood} on a 4-cell, 3-site GT10 example:
 * the Dirichlet-multinomial terms reproduce LikelihoodReadCountModel.logLiklihoodRC, and along MCMC-like
 * sequences of accepted and rejected read-count moves the likelihood equals one built from scratch,
 * with the terms recomputed only for epsilon, w1 and w2 and only the moved cell rebuilt for one size factor.
 * Genotype sampling through the likelihood matches the exact posterior marginals of the genotypes.
 */
public class ReadCountTreeLikelihoodCacheTest {

    private static final double DELTA = 1e-8;
    private static final int N = 10;
    private static final String[] TAXA = {"a", "b", "c", "d"};
    private static final String NEWICK = "((a:0.2,b:0.35):0.15,(c:0.1,d:0.4):0.25);";
    private static final String READS =
            "9:0:0:0,0:6:5:0,4:0:0:3\n"
            + "8:1:0:0,0:7:0:0,0:0:0:9\n"
            + "3:5:0:0,0:5:6:0,5:0:0:0\n"
            + "10:0:0:0,0:0:8:0,2:0:0:6";

    @BeforeAll
    public static void setUpClass() {
        addServices("version.xml");
    }

    private static RealScalarParam scalar(double v) {
        return new RealScalarParam(v, PositiveReal.INSTANCE);
    }

    private SiteModel siteModel() {
        double[] piVals = {0.3, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.1, 0.05, 0.1};
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", new SimplexParam(piVals), "estimate", false);
        RealVectorParam nucRates = new RealVectorParam(new double[]{1, 5, 1, 1, 5, 1}, PositiveReal.INSTANCE);
        nucRates.setInputValue("keys", "AC AG AT CG CT GT");
        nucRates.initAndValidate();
        GT10 subs = new GT10();
        subs.initByName("nucRates", nucRates, "frequencies", freqs);
        SiteModel sm = new SiteModel();
        sm.initByName("mutationRate", scalar(1.0), "gammaCategoryCount", 1, "substModel", subs);
        return sm;
    }

    private Alignment scaffold() {
        String[] seqs = {"ACG", "CGT", "GTA", "TAC"};
        Alignment data = new Alignment();
        for (int k = 0; k < TAXA.length; k++) data.setInputValue("sequence", new Sequence(TAXA[k], seqs[k]));
        data.setInputValue("dataType", "nucleotideDiploid10");
        data.initAndValidate();
        return data;
    }

    /** Read-count model; t and v small enough that coverage above 60 has negligible probability. */
    private LikelihoodReadCountModel readCountModel(Alignment scaffold) {
        ReadCount rc = new ReadCount();
        rc.initByName("value", READS, "taxaNames", String.join(" ", TAXA));
        LikelihoodReadCountModel rcm = new LikelihoodReadCountModel();
        rcm.initByName("alignment", scaffold, "readCount", rc,
                "epsilon", scalar(0.04), "delta", scalar(0.2), "t", scalar(3.0), "v", scalar(1.5),
                "s", new RealVectorParam(new double[]{1.0, 0.8, 1.2, 1.1}, PositiveReal.INSTANCE),
                "w1", scalar(5.0), "w2", scalar(3.0));
        return rcm;
    }

    /** Read-count model at given values {epsilon, delta, t, v, w1, w2} and size factors. */
    private LikelihoodReadCountModel readCountModel(Alignment scaffold, double[] p, double[] size) {
        ReadCount rc = new ReadCount();
        rc.initByName("value", READS, "taxaNames", String.join(" ", TAXA));
        LikelihoodReadCountModel rcm = new LikelihoodReadCountModel();
        rcm.initByName("alignment", scaffold, "readCount", rc,
                "epsilon", scalar(p[0]), "delta", scalar(p[1]), "t", scalar(p[2]), "v", scalar(p[3]),
                "s", new RealVectorParam(size.clone(), PositiveReal.INSTANCE),
                "w1", scalar(p[4]), "w2", scalar(p[5]));
        return rcm;
    }

    private ReadCountTreeLikelihood likelihoodAt(double[] p, double[] size, Object... extra) {
        Alignment scaffold = scaffold();
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", scaffold, "newick", NEWICK, "IsLabelledNewick", true);
        Object[] base = {"data", scaffold, "tree", tree, "siteModel", siteModel(),
                "readCountModel", readCountModel(scaffold, p, size)};
        Object[] args = new Object[base.length + extra.length];
        System.arraycopy(base, 0, args, 0, base.length);
        System.arraycopy(extra, 0, args, base.length, extra.length);
        ReadCountTreeLikelihood l = new ReadCountTreeLikelihood();
        l.initByName(args);
        return l;
    }

    private ReadCountTreeLikelihood likelihood(Object... extra) {
        Alignment scaffold = scaffold();
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", scaffold, "newick", NEWICK, "IsLabelledNewick", true);
        Object[] base = {"data", scaffold, "tree", tree, "siteModel", siteModel(),
                "readCountModel", readCountModel(scaffold)};
        Object[] args = new Object[base.length + extra.length];
        System.arraycopy(base, 0, args, 0, base.length);
        System.arraycopy(extra, 0, args, base.length, extra.length);
        ReadCountTreeLikelihood l = new ReadCountTreeLikelihood();
        l.initByName(args);
        return l;
    }

    private static final double[] SIZE = {1.0, 0.8, 1.2, 1.1};
    private static final String[] GT10 = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"};

    private static Frequencies germlineRoot(double het) {
        double[] f = new double[N];
        for (int g = 0; g < N; g++) f[g] = (g == 0 || g == 4 || g == 7 || g == 9) ? (1 - het) / 4 : het / 6;
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", new SimplexParam(f), "estimate", false);
        return freqs;
    }

    /**
     * A change to estimated root frequencies must trigger recalculation. In an MCMC run the State
     * calls CalculationNode.checkDirtiness() on the Frequencies node (which also refreshes its
     * cached values); that method is package-private, so the test invokes it by reflection.
     */
    @Test
    public void testRootFrequencyChangeIsSeen() throws Exception {
        Frequencies root = germlineRoot(0.01);
        ReadCountTreeLikelihood sel = likelihood("rootFrequencies", root);
        double logP0 = sel.calculateLogP();
        SimplexParam f = (SimplexParam) root.frequenciesInput.get();
        double[] v = new double[N];
        for (int g = 0; g < N; g++) v[g] = f.get(g);
        v[1] += 0.05;
        v[0] -= 0.05;
        sel.store();
        for (int g = 0; g < N; g++) f.set(g, v[g]);
        java.lang.reflect.Method check = beast.base.inference.CalculationNode.class.getDeclaredMethod("checkDirtiness");
        check.setAccessible(true);
        check.invoke(root);
        assertTrue(sel.requiresRecalculation());
        double logP1 = sel.calculateLogP();
        assertNotEquals(logP0, logP1, DELTA);
        ReadCountTreeLikelihood fresh = likelihood("rootFrequencies", root);
        assertEquals(fresh.calculateLogP(), logP1, DELTA);
    }

    // ------------------------------------------------ tip partials from cached Dirichlet terms

    @Test
    public void testTermsReproduceLogLikelihoodRC() {
        for (double[] p : new double[][]{{0.04, 0.2, 3.0, 1.5, 5.0, 3.0}, {0.003, 0.05, 6.0, 20.0, 0.7, 12.0},
                {1e-5, 1e-6, 0.5, 0.01, 50.0, 0.05}, {0.3, 0.999, 30.0, 200.0, 0.01, 300.0}}) {
            ReadCountTreeLikelihood l = likelihoodAt(p, SIZE);
            LikelihoodReadCountModel rcm = l.readCountModelInput.get();
            rcm.initialize();
            int nT = rcm.dirichletTermCount();
            double[] terms = new double[2 * nT + 1];
            double[] got = new double[N], tip = new double[N], scratch = new double[N];
            for (int c = 0; c < TAXA.length; c++) {
                for (int site = 0; site < 3; site++) {
                    int[] reads = rcm.getReadCount().getReadCounts(c, site);
                    int cov = rcm.getCoverage(c, site);
                    rcm.dirichletTerms(reads, cov, terms, 0);
                    rcm.logLikelihoodsFromTerms(terms, 0, cov, c, got);
                    rcm.linearTerms(terms, 0, terms, nT);
                    double norm = rcm.tipPartialsFromTerms(terms, 0, terms, nT, cov, c, tip, 0, scratch);
                    double max = Double.NEGATIVE_INFINITY;
                    for (int g = 0; g < N; g++) max = Math.max(max, rcm.logLiklihoodRC(g, reads, cov, c));
                    assertEquals(max, norm, 1e-10 * Math.max(1.0, Math.abs(max)), "normaliser, cell " + c + " site " + site);
                    for (int g = 0; g < N; g++) {
                        double ref = rcm.logLiklihoodRC(g, reads, cov, c);
                        assertEquals(ref, got[g], 0.0, "cell " + c + " site " + site + " g " + g);
                        // linear form: the same max-normalised partials
                        assertEquals(Math.exp(ref - max), tip[g], 1e-12, "partial, cell " + c + " site " + site + " g " + g);
                    }
                }
            }
        }
    }

    /**
     * MCMC-like sequences of read-count moves, accepted or rejected: after each, the likelihood equals
     * one built from scratch at the same values (so the cached terms and per-cell rebuilds are right).
     */
    @Test
    public void testCachedTipsFollowReadCountMoves() {
        Object[][] variants = {{}, {"rootFrequencies", germlineRoot(0.01)}};
        for (Object[] extra : variants) {
            double[] p = {0.04, 0.2, 3.0, 1.5, 5.0, 3.0};
            double[] size = SIZE.clone();
            ReadCountTreeLikelihood l = likelihoodAt(p, size, extra);
            LikelihoodReadCountModel rcm = l.readCountModelInput.get();
            RealScalarParam[] scal = {(RealScalarParam) rcm.epsilonInput.get(), (RealScalarParam) rcm.deltaInput.get(),
                    (RealScalarParam) rcm.tInput.get(), (RealScalarParam) rcm.vInput.get(),
                    (RealScalarParam) rcm.w1Input.get(), (RealScalarParam) rcm.w2Input.get()};
            RealVectorParam s = (RealVectorParam) rcm.sInput.get();
            beast.base.inference.State state = new beast.base.inference.State();
            state.initByName("stateNode", java.util.List.of(scal[0], scal[1], scal[2], scal[3], scal[4], scal[5], s));
            state.initialise();
            beast.base.inference.CompoundDistribution posterior = new beast.base.inference.CompoundDistribution();
            posterior.initByName("distribution", l);
            state.setPosterior(posterior);
            state.robustlyCalcPosterior(posterior);
            int step = 0;
            // move: parameter index (0-5 scalars, 6 one size factor, 7 all size factors), factor, accept
            Object[][] moves = {{0, 1.5, false}, {2, 1.2, true}, {6, 1.3, false}, {6, 0.8, true}, {0, 0.7, true},
                    {7, 1.1, false}, {1, 1.4, true}, {4, 0.9, false}, {6, 1.25, true}, {3, 1.3, true}, {5, 0.8, false},
                    {0, 1.2, false}, {6, 0.9, true}};
            int cell = 0;
            for (Object[] mv : moves) {
                int which = (Integer) mv[0];
                double f = (Double) mv[1];
                boolean accept = (Boolean) mv[2];
                double before = l.calculateLogP();
                double[] p0 = p.clone(), size0 = size.clone();
                // MCMC's order: store the state, propose, then store and check the calculation nodes
                state.store(++step);
                if (which < 6) {
                    p[which] *= f;
                    scal[which].set(p[which]);
                } else if (which == 6) {
                    cell = (cell + 1) % TAXA.length;
                    size[cell] *= f;
                    s.set(cell, size[cell]);
                } else {
                    for (int c = 0; c < TAXA.length; c++) { size[c] *= f; s.set(c, size[c]); }
                }
                state.storeCalculationNodes();
                state.checkCalculationNodesDirtiness();
                int refills = l.termRefills;
                double moved = l.calculateLogP();
                // only epsilon, w1, w2 recompute the Dirichlet terms; one size factor rebuilds one cell
                assertEquals(which == 0 || which == 4 || which == 5 ? refills + 1 : refills, l.termRefills,
                        "term refills after move " + which);
                assertEquals(which == 6 ? 1 : TAXA.length, l.leavesRebuilt, "leaves rebuilt after move " + which);
                assertEquals(likelihoodAt(p, size, extra).calculateLogP(), moved, 1e-9, "after move " + which + " (step "
                        + step + ", before " + before + ", " + java.util.Arrays.toString(extra) + ")");
                if (accept) {
                    state.acceptCalculationNodes();
                } else {
                    state.restore();
                    state.restoreCalculationNodes();
                    p = p0;
                    size = size0;
                    assertEquals(before, l.calculateLogP(), 0.0, "restored after move " + which);
                }
                state.setEverythingDirty(false);   // as MCMC does at the end of each step
                assertEquals(likelihoodAt(p, size, extra).calculateLogP(), l.calculateLogP(), 1e-9,
                        "state after move " + which + (accept ? " accepted" : " rejected"));
            }
        }
    }

    // ------------------------------------------------ genotype sampling

    /**
     * sampleGenotypes draws from the posterior of the genotypes given the reads: the sampled marginals
     * of every leaf and of the root match the exact ones (summing over the root and both internal
     * nodes of the 4-cell tree).
     */
    @Test
    public void testSampleGenotypesMatchesExactMarginals() {
        ReadCountTreeLikelihood l = likelihood("rootFrequencies", germlineRoot(0.01));
        l.calculateLogP();
        TreeParser tree = (TreeParser) l.treeInput.get();
        GT10 subs = (GT10) ((SiteModel) l.siteModelInput.get()).substModelInput.get();
        LikelihoodReadCountModel rcm = l.readCountModelInput.get();
        rcm.initialize();
        Node root = tree.getRoot(), x = root.getLeft(), y = root.getRight();
        double[][] P = new double[tree.getNodeCount()][N * N];
        for (Node node : tree.getNodesAsArray()) {
            if (!node.isRoot()) subs.getTransitionProbabilities(node, node.getLength(), 0.0, 1.0, P[node.getNr()]);
        }
        double[] pi = germlineRoot(0.01).getFreqs();
        int nSamples = 20000;
        int[][] states = new int[3][tree.getNodeCount()];
        java.util.Random random = new java.util.Random(3);
        int[][][] counts = new int[3][tree.getNodeCount()][N];
        for (int i = 0; i < nSamples; i++) {
            l.sampleGenotypes(random::nextDouble, states);
            for (int k = 0; k < 3; k++) for (int nr = 0; nr < tree.getNodeCount(); nr++) counts[k][nr][states[k][nr]]++;
        }
        for (int k = 0; k < 3; k++) {
            double[] rootDist = pi;
            // tip likelihoods
            double[][] tip = new double[tree.getNodeCount()][N];
            for (Node leaf : tree.getExternalNodes()) {
                int c = java.util.Arrays.asList(TAXA).indexOf(leaf.getID());
                int[] reads = rcm.getReadCount().getReadCounts(c, k);
                for (int g = 0; g < N; g++) tip[leaf.getNr()][g] = Math.exp(rcm.logLiklihoodRC(g, reads, rcm.getCoverage(c, k), c));
            }
            double[][] exact = new double[tree.getNodeCount()][N];
            double total = 0.0;
            for (int r = 0; r < N; r++) for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) {
                double w = rootDist[r] * P[x.getNr()][r * N + a] * P[y.getNr()][r * N + b];
                // leaf terms below x and y, each leaf's own weight vector, and their sums
                Node[] leaves = {x.getLeft(), x.getRight(), y.getLeft(), y.getRight()};
                int[] parent = {a, a, b, b};
                double[][] lw = new double[4][N];
                double[] ls = new double[4];
                for (int j = 0; j < 4; j++) {
                    for (int g = 0; g < N; g++) {
                        lw[j][g] = P[leaves[j].getNr()][parent[j] * N + g] * tip[leaves[j].getNr()][g];
                        ls[j] += lw[j][g];
                    }
                }
                double joint = w * ls[0] * ls[1] * ls[2] * ls[3];
                total += joint;
                exact[root.getNr()][r] += joint;
                for (int j = 0; j < 4; j++) {
                    for (int g = 0; g < N; g++) exact[leaves[j].getNr()][g] += joint / ls[j] * lw[j][g];
                }
            }
            for (int nr : new int[]{root.getNr(), x.getLeft().getNr(), x.getRight().getNr(), y.getLeft().getNr(),
                    y.getRight().getNr()}) {
                for (int g = 0; g < N; g++) {
                    double p = exact[nr][g] / total, f = counts[k][nr][g] / (double) nSamples;
                    double tol = 5 * Math.sqrt(p * (1 - p) / nSamples) + 1e-3;
                    assertEquals(p, f, tol, "site " + k + " node " + nr + " genotype " + GT10[g]);
                }
            }
        }
    }

    /** The logger finds the likelihood that uses its read-count model, and logs one genotype per site. */
    @Test
    public void testGenotypeLoggerSamplesThroughTheLikelihood() {
        ReadCountTreeLikelihood l = likelihood("rootFrequencies", germlineRoot(0.01));
        l.calculateLogP();
        LikelihoodReadCountModel rcm = l.readCountModelInput.get();
        phylonco.beast.evolution.logger.SampledGenotypeLogger logger = new phylonco.beast.evolution.logger.SampledGenotypeLogger();
        logger.initByName("tree", l.treeInput.get(), "siteModel", l.siteModelInput.get(), "readCountModel", rcm,
                "readCount", rcm.getReadCount());
        assertTrue(logger.likelihood() == l, "logger samples through the ReadCountTreeLikelihood");
        java.io.ByteArrayOutputStream bytes = new java.io.ByteArrayOutputStream();
        java.io.PrintStream out = new java.io.PrintStream(bytes);
        logger.init(out);
        out.println();
        logger.log(0, out);
        out.flush();
        String[] lines = bytes.toString().trim().split("\n");
        String[] header = lines[0].trim().split("\t"), row = lines[1].trim().split("\t");
        assertEquals(TAXA.length, header.length);
        assertEquals(TAXA.length, row.length);
        for (String g : row) assertEquals(3, g.length(), "one genotype per site: " + g);
    }
}
