package phylonco.beast.evolution.operator;

import beast.base.evolution.operator.kernel.BactrianSubtreeSlide;
import beast.base.inference.Operator;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.evolution.operator.Uniform;
import beast.base.evolution.speciation.YuleModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.TreeStatLogger;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * Tests for {@link BactrianSubtreeScale}.
 *
 * <p>Two checks:</p>
 * <ol>
 *   <li><b>Mechanics</b> — on an ultrametric tree the operator never proposes a node below
 *       its carried subtree (no negative branches) and never flat-rejects from that bound.</li>
 *   <li><b>Detailed balance</b> — sampling a Yule prior with this operator reproduces the same
 *       tree-height distribution as the trusted {@link BactrianSubtreeSlide}; a wrong Hastings
 *       ratio would bias the height.</li>
 * </ol>
 */
public class BactrianSubtreeScaleTest {

    private static final String NEWICK =
            "((A:2,B:2):1,(C:2,(D:1,(E:0.5,F:0.5):0.5):1):1):0;";

    private Tree newTree() {
        Tree tree = new TreeParser();
        tree.initByName(
                "newick", NEWICK,
                "IsLabelledNewick", true,
                "adjustTipHeights", false);
        return tree;
    }

    /** Every non-root node must sit strictly below its parent (no zero/negative branches). */
    private void assertValid(Tree tree) {
        for (Node node : tree.getNodesAsArray()) {
            if (!node.isRoot()) {
                assertTrue("node " + node.getNr() + " height " + node.getHeight()
                                + " not below parent " + node.getParent().getHeight(),
                        node.getHeight() < node.getParent().getHeight() + 1e-12);
            }
            assertTrue("negative height", node.getHeight() >= 0.0);
        }
    }

    @Test
    public void testFloorNeverViolated() {
        Randomizer.setSeed(42);
        Tree tree = newTree();

        BactrianSubtreeScale op = new BactrianSubtreeScale();
        op.initByName("tree", tree, "size", 0.5, "weight", 1.0);

        int infinite = 0;
        final int N = 20000;
        for (int k = 0; k < N; k++) {
            double logHR = op.proposal();
            if (logHR == Double.NEGATIVE_INFINITY) {
                infinite++;
                // tree must be untouched on a rejected-before-mutation proposal
                assertValid(tree);
            } else {
                assertTrue("non-finite logHR " + logHR, !Double.isNaN(logHR) && !Double.isInfinite(logHR));
                assertValid(tree);
            }
        }
        // on an ultrametric tree the carried-subtree floor and the destination check can never fire
        assertTrue("unexpected flat rejects on ultrametric tree: " + infinite, infinite == 0);
    }

    /** Mean tree height over the post-burnin samples of a 2-column (Sample, height) trace file. */
    private double meanHeight(Path log) throws IOException {
        List<String> lines = Files.readAllLines(log);
        List<Double> heights = new ArrayList<>();
        for (String line : lines) {
            if (line.startsWith("Sample") || line.startsWith("#") || line.trim().isEmpty()) continue;
            String[] f = line.split("\\s+");
            if (f.length < 2) continue;
            try {
                heights.add(Double.parseDouble(f[1]));
            } catch (NumberFormatException ignore) {
                // header / non-numeric
            }
        }
        int burnin = heights.size() / 10;
        double sum = 0;
        int n = 0;
        for (int k = burnin; k < heights.size(); k++) {
            sum += heights.get(k);
            n++;
        }
        return sum / n;
    }

    private double sampleYuleMeanHeight(Operator treeMove, Tree tree, long seed) throws Exception {
        Randomizer.setSeed(seed);

        RealParameter birthRate = new RealParameter("1.0");

        YuleModel yule = new YuleModel();
        yule.initByName("tree", tree, "birthDiffRate", birthRate);

        State state = new State();
        state.initByName("stateNode", tree);

        CompoundDistribution posterior = new CompoundDistribution();
        posterior.initByName("distribution", (Distribution) yule);

        Uniform uniform = new Uniform();
        uniform.initByName("tree", tree, "weight", 3.0);

        ScaleOperator scale = new ScaleOperator();
        scale.initByName("tree", tree, "scaleFactor", 0.75, "weight", 1.0);

        List<Operator> operators = new ArrayList<>();
        operators.add(uniform);
        operators.add(scale);
        operators.add(treeMove);

        Path logFile = Files.createTempFile("subtreeScaleTest", ".log");
        Files.delete(logFile); // let BEAST create it fresh

        TreeStatLogger heightLogger = new TreeStatLogger();
        heightLogger.initByName("tree", tree, "logLength", false);

        Logger logger = new Logger();
        logger.initByName("logEvery", 2000, "fileName", logFile.toString(), "log", heightLogger);

        MCMC mcmc = new MCMC();
        mcmc.initByName(
                "chainLength", 6000000L,
                "state", state,
                "distribution", posterior,
                "operator", operators,
                "logger", logger);

        mcmc.run();

        return meanHeight(logFile);
    }

    @Test
    public void testYuleHeightMatchesSubtreeSlide() throws Exception {
        // reference: trusted additive subtree slide
        Tree refTree = newTree();
        BactrianSubtreeSlide slide = new BactrianSubtreeSlide();
        slide.initByName("tree", refTree, "weight", 5.0);
        double meanSlide = sampleYuleMeanHeight(slide, refTree, 123);

        // operator under test
        Tree testTree = newTree();
        BactrianSubtreeScale scaleMove = new BactrianSubtreeScale();
        scaleMove.initByName("tree", testTree, "size", 0.5, "weight", 5.0);
        double meanScale = sampleYuleMeanHeight(scaleMove, testTree, 321);

        double relDiff = Math.abs(meanScale - meanSlide) / meanSlide;
        assertTrue("Yule mean tree height disagrees: slide=" + meanSlide
                + " scale=" + meanScale + " relDiff=" + relDiff, relDiff < 0.05);
    }
}
