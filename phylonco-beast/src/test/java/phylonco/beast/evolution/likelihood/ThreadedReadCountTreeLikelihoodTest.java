package phylonco.beast.evolution.likelihood;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.State;
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
import java.util.List;
import java.util.Random;

import static beast.pkgmgmt.BEASTClassLoader.addServices;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * {@link ThreadedReadCountTreeLikelihood} gives exactly the log likelihood of
 * {@link ReadCountTreeLikelihood}, whatever the number of threads, at the start and along MCMC-like
 * sequences of accepted and rejected tree and read-count moves.
 */
public class ThreadedReadCountTreeLikelihoodTest {

    private static final String[] TAXA = {"a", "b", "c", "d", "e", "f", "g"};
    private static final String NEWICK =
            "(((a:0.1,b:0.25):0.2,(c:0.15,d:0.05):0.3):0.1,((e:0.2,f:0.35):0.1,g:0.45):0.15);";
    private static final int L = 60;
    private static final String READS = randomReads(new Random(17));

    @BeforeAll
    public static void setUpClass() {
        addServices("version.xml");
    }

    private static String randomReads(Random r) {
        StringBuilder sb = new StringBuilder();
        for (int t = 0; t < TAXA.length; t++) {
            if (t > 0) sb.append("\n");
            for (int j = 0; j < L; j++) {
                if (j > 0) sb.append(",");
                int[] c = new int[4];
                int ref = r.nextInt(4), cov = r.nextInt(25);
                for (int k = 0; k < cov; k++) c[r.nextDouble() < 0.7 ? ref : r.nextInt(4)]++;
                sb.append(c[0]).append(":").append(c[1]).append(":").append(c[2]).append(":").append(c[3]);
            }
        }
        return sb.toString();
    }

    private static RealScalarParam scalar(double v) {
        return new RealScalarParam(v, PositiveReal.INSTANCE);
    }

    /** One model graph: its likelihood, tree, read-count parameters and a State over them. */
    private static final class Model {
        ReadCountTreeLikelihood l;
        TreeParser tree;
        RealScalarParam[] scal;
        RealVectorParam s;
        State state;
    }

    private Model build(int threads, Object... extra) {
        return build(threads, false, extra);
    }

    /** threads 0: ReadCountTreeLikelihood; above 0: the threaded subclass; beagle: BeagleReadCountTreeLikelihood. */
    private Model build(int threads, boolean beagle, Object... extra) {
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

        Alignment scaffold = new Alignment();
        // distinct columns (the likelihood needs one pattern per site; the genotype values are ignored):
        // column j has digit t of j, base 10, as taxon t's genotype
        String codes = "ACGTMRWSYK";
        for (int t = 0; t < TAXA.length; t++) {
            StringBuilder seq = new StringBuilder();
            for (int j = 0; j < L; j++) seq.append(codes.charAt((int) (j / Math.pow(10, t)) % 10));
            scaffold.setInputValue("sequence", new Sequence(TAXA[t], seq.toString()));
        }
        scaffold.setInputValue("dataType", "nucleotideDiploid10");
        scaffold.initAndValidate();
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", scaffold, "newick", NEWICK, "IsLabelledNewick", true);

        ReadCount rc = new ReadCount();
        rc.initByName("value", READS, "taxaNames", String.join(" ", TAXA));
        Model m = new Model();
        m.scal = new RealScalarParam[]{scalar(0.03), scalar(0.15), scalar(4.0), scalar(2.0), scalar(5.0), scalar(3.0)};
        double[] size = new double[TAXA.length];
        for (int i = 0; i < size.length; i++) size[i] = 0.8 + 0.07 * i;
        m.s = new RealVectorParam(size, PositiveReal.INSTANCE);
        LikelihoodReadCountModel rcm = new LikelihoodReadCountModel();
        rcm.initByName("alignment", scaffold, "readCount", rc, "epsilon", m.scal[0], "delta", m.scal[1],
                "t", m.scal[2], "v", m.scal[3], "s", m.s, "w1", m.scal[4], "w2", m.scal[5]);

        Object[] base = {"data", scaffold, "tree", tree, "siteModel", sm, "readCountModel", rcm};
        Object[] args = Arrays.copyOf(base, base.length + extra.length + (threads > 0 ? 2 : 0));
        System.arraycopy(extra, 0, args, base.length, extra.length);
        if (threads > 0) {
            args[args.length - 2] = "threads";
            args[args.length - 1] = threads;
        }
        m.l = beagle ? new BeagleReadCountTreeLikelihood()
                : threads > 0 ? new ThreadedReadCountTreeLikelihood() : new ReadCountTreeLikelihood();
        m.l.initByName(args);
        m.tree = tree;

        m.state = new State();
        m.state.initByName("stateNode", List.of(tree, m.scal[0], m.scal[1], m.scal[2], m.scal[3], m.scal[4],
                m.scal[5], m.s));
        m.state.initialise();
        CompoundDistribution posterior = new CompoundDistribution();
        posterior.initByName("distribution", m.l);
        m.state.setPosterior(posterior);
        m.state.robustlyCalcPosterior(posterior);
        return m;
    }

    private static final Object[][] VARIANTS = {
            {},
            {"rootFrequencies", null}};

    private static Object[] variant(int v) {
        Object[] x = VARIANTS[v].clone();
        for (int i = 0; i < x.length; i++) {
            if ("rootFrequencies".equals(x[i])) {
                double[] f = new double[10];
                for (int g = 0; g < 10; g++) f[g] = (g == 0 || g == 4 || g == 7 || g == 9) ? 0.2475 : 0.01 / 6;
                Frequencies root = new Frequencies();
                root.initByName("frequencies", new SimplexParam(f), "estimate", false);
                x[i + 1] = root;
            }
        }
        return x;
    }

    @Test
    public void testSameLogLikelihoodForAnyThreadCount() {
        for (int v = 0; v < VARIANTS.length; v++) {
            double expected = build(0, variant(v)).l.calculateLogP();
            for (int threads : new int[]{1, 2, 3, 8}) {
                assertEquals(expected, build(threads, variant(v)).l.calculateLogP(), 0.0,
                        "variant " + v + ", " + threads + " threads");
            }
        }
    }

    private static boolean beagleAvailable() {
        try {
            return !beagle.BeagleFactory.getResourceDetails().isEmpty();
        } catch (Throwable e) {
            return false;
        }
    }

    /**
     * BeagleReadCountTreeLikelihood gives ReadCountTreeLikelihood's log likelihood (to rounding: BEAGLE
     * rescales every node) for every variant; skipped without the BEAGLE native library.
     */
    @Test
    public void testBeagleMatchesReadCountTreeLikelihood() {
        org.junit.jupiter.api.Assumptions.assumeTrue(beagleAvailable(), "BEAGLE native library not available");
        for (int v = 0; v < VARIANTS.length; v++) {
            double java = build(0, variant(v)).l.calculateLogP();
            for (int threads : new int[]{0, 2}) {
                double b = build(threads, true, variant(v)).l.calculateLogP();
                assertEquals(java, b, 1e-10 * Math.abs(java), "variant " + v + ", " + threads + " threads");
            }
        }
    }

    /** As testThreadedFollowsMovesExactly, BEAGLE against ReadCountTreeLikelihood, to rounding. */
    @Test
    public void testBeagleFollowsMoves() {
        org.junit.jupiter.api.Assumptions.assumeTrue(beagleAvailable(), "BEAGLE native library not available");
        followMoves(true, 1e-10);
    }

    /** Apply one move to a model: 0-5 a read-count scalar, 6 one size factor, 7 every size factor, 8 a node height. */
    private static void move(Model m, int which, double f, int cell, int nodeNr) {
        if (which < 6) {
            m.scal[which].set(m.scal[which].get() * f);
        } else if (which == 6) {
            m.s.set(cell, m.s.get(cell) * f);
        } else if (which == 7) {
            for (int c = 0; c < TAXA.length; c++) m.s.set(c, m.s.get(c) * f);
        } else {
            Node n = m.tree.getNode(nodeNr);
            double lo = Math.max(n.getLeft().getHeight(), n.getRight().getHeight());
            double hi = n.getParent().getHeight();
            n.setHeight(lo + f * (hi - lo));
        }
    }

    @Test
    public void testThreadedFollowsMovesExactly() {
        followMoves(false, 0.0);
    }

    private void followMoves(boolean beagle, double relTol) {
        // move, factor, accept
        Object[][] moves = {{8, 0.3, true}, {0, 1.4, false}, {6, 1.2, true}, {8, 0.7, false}, {2, 0.9, true},
                {7, 1.1, false}, {4, 0.8, true}, {8, 0.55, true}, {1, 1.3, false}, {6, 0.85, false}, {3, 1.2, true},
                {0, 0.8, true}, {8, 0.2, false}, {5, 1.25, true}, {6, 1.1, true}};
        for (int v = 0; v < VARIANTS.length; v++) {
            for (int threads : new int[]{2, 3}) {
                Model ref = build(0, variant(v)), thr = build(threads, beagle, variant(v));
                int step = 0, cell = 0, internal = 0;
                List<Node> inner = new java.util.ArrayList<>();
                for (Node n : ref.tree.getInternalNodes()) if (!n.isRoot()) inner.add(n);
                for (Object[] mv : moves) {
                    int which = (Integer) mv[0];
                    double f = (Double) mv[1];
                    boolean accept = (Boolean) mv[2];
                    cell = (cell + 3) % TAXA.length;
                    internal = (internal + 1) % inner.size();
                    int nodeNr = inner.get(internal).getNr();
                    step++;
                    for (Model m : new Model[]{ref, thr}) {
                        m.state.store(step);
                        move(m, which, f, cell, nodeNr);
                        m.state.storeCalculationNodes();
                        m.state.checkCalculationNodesDirtiness();
                    }
                    double a = ref.l.calculateLogP(), b = thr.l.calculateLogP();
                    assertTrue(Double.isFinite(a), "finite log likelihood");
                    assertEquals(a, b, relTol * Math.abs(a), "variant " + v + ", " + threads + " threads, move " + step);
                    for (Model m : new Model[]{ref, thr}) {
                        if (accept) {
                            m.state.acceptCalculationNodes();
                        } else {
                            m.state.restore();
                            m.state.restoreCalculationNodes();
                        }
                        m.state.setEverythingDirty(false);
                    }
                    double ra = ref.l.calculateLogP();
                    assertEquals(ra, thr.l.calculateLogP(), relTol * Math.abs(ra),
                            "variant " + v + ", " + threads + " threads, after move " + step);
                }
            }
        }
    }
}
