package phylonco.beast.evolution.operator;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

/**
 * Subtree-slide move that proposes the moved node's height by scaling, in log
 * space, the branch between the node and the subtree it carries.
 *
 * <p>When a node {@code p} is slid, it carries the subtree rooted at one of its
 * children ({@code i}); that subtree's heights do not change, so {@code i}'s
 * height is a hard floor for {@code p}. {@link beast.base.evolution.operator.kernel.BactrianSubtreeSlide}
 * proposes {@code p}'s height additively, so a downward step routinely overshoots
 * this floor and is flat-rejected ({@code NEGATIVE_INFINITY}); with a large step
 * size that flat-rejection dominates and the operator cannot be auto-tuned into a
 * useful regime.</p>
 *
 * <p>This operator instead proposes
 * <pre>
 *   gap'    = gap * exp(size * kernel)        // kernel is a mean-0 Bactrian draw
 *   height' = floor + gap'
 * </pre>
 * where {@code gap = p.height - floor > 0}. Because {@code gap' > 0} for any draw,
 * the move can never put {@code p} below the carried subtree, so that source of
 * rejection is removed by construction. The proposal is also scale-invariant with
 * respect to branch length, so a single {@code size} works across branches of very
 * different lengths and the operator auto-optimises cleanly.</p>
 *
 * <h3>Hastings ratio</h3>
 * <p>The proposal factorises into the continuous height move and the discrete
 * reattachment, so {@code logHR = log(scale) + combinatorial}, where
 * {@code log(scale) = log(gap'/gap)} is the standard scale-operator Jacobian and the
 * combinatorial term ({@code -log(#sources)} up, {@code +log(#destinations)} down,
 * {@code 0} for a pure node-age move) is identical to {@code BactrianSubtreeSlide}.</p>
 *
 * <p>Rejections from the destination side remain (e.g. in non-ultrametric trees the
 * proposed height can fall below every leaf of the reattachment subtree, leaving no
 * straddling edge); these are inherent to subtree slide and are not affected here.</p>
 */
@Description("Subtree-slide move that scales the branch above the carried subtree in log space, "
        + "so the moved node can never be proposed below the carried subtree (removing that flat rejection) "
        + "and the step is scale-invariant w.r.t. branch length, making it easier to auto-optimise than "
        + "an additive subtree slide.")
public class BactrianSubtreeScale extends TreeOperator {

    final public Input<Double> sizeInput = new Input<>("size",
            "scale of the log-space random walk on the gap above the carried subtree, default 0.3", 0.3);
    final public Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution",
            "provides sample distribution for proposals", KernelDistribution.newDefaultKernelDistribution());
    final public Input<Boolean> optimiseInput = new Input<>("optimise",
            "flag to indicate that the size is automatically changed to achieve a good acceptance rate (default true)", true);

    private double size;
    private KernelDistribution kernelDistribution;

    @Override
    public void initAndValidate() {
        size = sizeInput.get();
        kernelDistribution = kernelDistributionInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = (Tree) InputUtil.get(treeInput, this);

        // 1. choose a random node avoiding the root; p is the node being slid.
        final int nodeCount = tree.getNodeCount();
        Node i;
        do {
            i = tree.getNode(Randomizer.nextInt(nodeCount));
        } while (i.isRoot());

        final Node p = i.getParent();
        final Node CiP = getOtherChild(p, i);   // sibling that p slides past / into
        final Node PiP = p.getParent();

        // 2. carried subtree (rooted at i) sets the floor; scale the gap above it.
        final double floor = i.getHeight();
        final double oldHeight = p.getHeight();
        final double oldGap = oldHeight - floor;
        if (oldGap <= 0.0) {
            // zero-length branch (e.g. sampled ancestor): nothing to scale
            return Double.NEGATIVE_INFINITY;
        }

        final double logScale = kernelDistribution.getRandomDelta(0, Double.NaN, size); // size * mean-0 draw
        final double scale = Math.exp(logScale);
        final double newHeight = floor + oldGap * scale;

        // Jacobian of the scale move on the gap.
        double logq = logScale;

        // 3. if the move is up
        if (newHeight > oldHeight) {

            // 3.1 if the topology will change
            if (PiP != null && PiP.getHeight() < newHeight) {
                // find new parent
                Node newParent = PiP;
                Node newChild = p;
                while (newParent.getHeight() < newHeight) {
                    newChild = newParent;
                    newParent = newParent.getParent();
                    if (newParent == null) break;
                }

                // 3.1.1 if creating a new root
                if (newChild.isRoot()) {
                    replace(p, CiP, newChild);
                    replace(PiP, p, CiP);

                    p.setParent(null);
                    tree.setRoot(p);
                }
                // 3.1.2 no new root
                else {
                    replace(p, CiP, newChild);
                    replace(PiP, p, CiP);
                    replace(newParent, newChild, p);
                }

                p.setHeight(newHeight);

                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(newChild, oldHeight, null);
                logq -= Math.log(possibleSources);

            } else {
                // just change the node height
                p.setHeight(newHeight);
            }
        }
        // 4. if we are sliding the subtree down
        else {

            // 4.0 is it a valid move? By construction newHeight > floor = i.getHeight(),
            // so the carried-subtree bound can never be violated. Guard defensively.
            if (i.getHeight() > newHeight) {
                return Double.NEGATIVE_INFINITY;
            }

            // 4.1 will the move change the topology
            if (CiP.getHeight() > newHeight) {

                final List<Node> newChildren = new ArrayList<>();
                final int possibleDestinations = intersectingEdges(CiP, newHeight, newChildren);

                // if no valid destinations then return a failure (e.g. non-ultrametric leaf ages)
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                final int childIndex = Randomizer.nextInt(newChildren.size());
                final Node newChild = newChildren.get(childIndex);
                final Node newParent = newChild.getParent();

                // 4.1.1 if p was root
                if (p.isRoot()) {
                    replace(p, CiP, newChild);
                    replace(newParent, newChild, p);

                    CiP.setParent(null);
                    tree.setRoot(CiP);

                } else {
                    replace(p, CiP, newChild);
                    replace(PiP, p, CiP);
                    replace(newParent, newChild, p);
                }

                p.setHeight(newHeight);

                logq += Math.log(possibleDestinations);
            } else {
                p.setHeight(newHeight);
            }
        }
        return logq;
    }

    protected int intersectingEdges(Node node, double height, List<Node> directChildren) {
        final Node parent = node.getParent();

        if (parent.getHeight() < height) return 0;

        if (node.getHeight() < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }

        if (node.isLeaf()) {
            return 0;
        } else {
            return intersectingEdges(node.getLeft(), height, directChildren) +
                    intersectingEdges(node.getRight(), height, directChildren);
        }
    }

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            // size is a Bactrian scale: larger size means larger jumps. Tune multiplicatively
            // so acceptance too high (delta > 0) increases size, too low decreases it.
            double delta = calcDelta(logAlpha);
            delta += Math.log(size);
            size = Math.exp(delta);
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return size;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        size = value;
    }

    @Override
    public double getTargetAcceptanceProbability() {
        return 0.3;
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newSize = size * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newSize);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newSize);
        } else return "";
    }
}
