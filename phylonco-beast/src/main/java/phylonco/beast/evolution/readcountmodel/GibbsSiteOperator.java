package phylonco.beast.evolution.readcountmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import mutablealignment.MutableAlignment;
import phylonco.beast.evolution.datatype.ReadCount;

import java.util.Arrays;

/**
 * Column-wise (site-wise) Gibbs sampler for leaf genotypes.
 *
 * Samples all leaf genotypes at a single site simultaneously from
 * P(leaf_states_at_site | tree, read_counts, parameters).
 *
 * Uses Felsenstein's pruning algorithm with read count likelihoods at leaves,
 * then samples from root down using the computed partials.
 * Internal nodes are sampled transiently to condition the leaf sampling properly.
 *
 * <p>Example XML usage (sample one random site per proposal):</p>
 * <pre>
 * &lt;operator id="GibbsSiteOperator" spec="phylonco.beast.evolution.readcountmodel.GibbsSiteOperator"
 *           weight="1.0"
 *           mutableAlignment="@A"
 *           tree="@psi"
 *           siteModel="@siteModel"
 *           readCountModel="@readCountLikelihood"
 *           readCount="@readCounts"/&gt;
 * </pre>
 *
 * <p>Example XML usage (sample all sites per proposal):</p>
 * <pre>
 * &lt;operator id="GibbsAlignmentSiteOperator" spec="phylonco.beast.evolution.readcountmodel.GibbsSiteOperator"
 *           weight="1.0"
 *           sampleAllSites="true"
 *           mutableAlignment="@A"
 *           tree="@psi"
 *           siteModel="@siteModel"
 *           readCountModel="@readCountLikelihood"
 *           readCount="@readCounts"/&gt;
 * </pre>
 *
 * <p>Where:</p>
 * <ul>
 *   <li>mutableAlignment: reference to a MutableAlignment containing genotype states</li>
 *   <li>tree: reference to the phylogenetic tree</li>
 *   <li>siteModel: reference to the SiteModel (contains substitution model)</li>
 *   <li>readCountModel: reference to LikelihoodReadCountModel</li>
 *   <li>readCount: reference to ReadCount data</li>
 *   <li>sampleAllSites: (optional) if true, sample all sites per proposal; default false</li>
 *   <li>branchRateModel: (optional) reference to BranchRateModel for relaxed clocks</li>
 * </ul>
 */
@Description("Gibbs sampler that samples all leaf genotypes at a single site, " +
        "using read count likelihoods and tree structure")
public class GibbsSiteOperator extends Operator {

    public Input<MutableAlignment> mutableAlignmentInput = new Input<>(
            "mutableAlignment",
            "mutable alignment containing genotype states",
            Input.Validate.REQUIRED);

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "phylogenetic tree",
            Input.Validate.REQUIRED);

    public Input<SiteModel.Base> siteModelInput = new Input<>(
            "siteModel",
            "site model containing substitution model",
            Input.Validate.REQUIRED);

    public Input<BranchRateModel.Base> branchRateModelInput = new Input<>(
            "branchRateModel",
            "branch rate model for relaxed clock (optional)",
            Input.Validate.OPTIONAL);

    public Input<LikelihoodReadCountModel> readCountModelInput = new Input<>(
            "readCountModel",
            "read count likelihood model",
            Input.Validate.REQUIRED);

    public Input<ReadCount> readCountInput = new Input<>(
            "readCount",
            "read count data",
            Input.Validate.REQUIRED);

    public Input<Boolean> sampleAllSitesInput = new Input<>(
            "sampleAllSites",
            "if true, sample all sites in one proposal; if false (default), sample one random site",
            false);

    // Cached references
    private MutableAlignment alignment;
    private Tree tree;
    private SubstitutionModel substitutionModel;
    private SiteModel.Base siteModel;
    private BranchRateModel.Base branchRateModel;
    private LikelihoodReadCountModel readCountModel;
    private ReadCount readCount;
    private boolean sampleAllSites;

    // Dimensions
    private int numNodes;
    private int numStates;
    private int numSites;
    private int numTaxa;

    // Pre-allocated arrays for efficiency
    private double[][] partials;        // [numNodes][numStates] - partial likelihoods
    private int[] sampledStates;        // [numNodes] - sampled state for each node
    private int[] leafStates;           // [numTaxa] - leaf states for alignment update
    private double[] transitionMatrix;  // [numStates * numStates] - reusable matrix
    private double[] logProbs;          // [numStates] - reusable probability array
    private int[][] coverages;          // [numTaxa][numSites] - precomputed coverages
    private double[] rootFrequencies;   // [numStates] - equilibrium frequencies

    // Mapping from node number to taxon index
    private int[] nodeNrToTaxonIndex;

    @Override
    public void initAndValidate() {
        alignment = mutableAlignmentInput.get();
        tree = treeInput.get();
        siteModel = siteModelInput.get();
        substitutionModel = siteModel.getSubstitutionModel();
        branchRateModel = branchRateModelInput.get();
        readCountModel = readCountModelInput.get();
        readCount = readCountInput.get();
        sampleAllSites = sampleAllSitesInput.get();

        numStates = alignment.getDataType().getStateCount();
        numSites = alignment.getSiteCount();
        numTaxa = alignment.getTaxonCount();
        numNodes = tree.getNodeCount();

        // Pre-allocate arrays
        partials = new double[numNodes][numStates];
        sampledStates = new int[numNodes];
        leafStates = new int[numTaxa];
        transitionMatrix = new double[numStates * numStates];
        logProbs = new double[numStates];

        // Precompute coverages for all taxa and sites
        coverages = new int[numTaxa][numSites];
        for (int i = 0; i < numTaxa; i++) {
            for (int j = 0; j < numSites; j++) {
                int[] counts = readCount.getReadCounts(i, j);
                for (int k = 0; k < counts.length; k++) {
                    coverages[i][j] += counts[k];
                }
            }
        }

        // Build mapping from node numbers to taxon indices
        nodeNrToTaxonIndex = new int[numNodes];
        Arrays.fill(nodeNrToTaxonIndex, -1);
        for (Node leaf : tree.getExternalNodes()) {
            String taxonName = leaf.getID();
            int taxonIndex = alignment.getTaxonIndex(taxonName);
            nodeNrToTaxonIndex[leaf.getNr()] = taxonIndex;
        }
    }

    @Override
    public double proposal() {
        if (sampleAllSites) {
            // Sample all sites
            for (int siteIndex = 0; siteIndex < numSites; siteIndex++) {
                sampleSite(siteIndex);
            }
        } else {
            // Sample one random site
            int siteIndex = Randomizer.nextInt(numSites);
            sampleSite(siteIndex);
        }

        // Gibbs operator: always accept
        return Double.POSITIVE_INFINITY;
    }

    /**
     * Sample all leaf genotypes at a single site.
     */
    private void sampleSite(int siteIndex) {
        // 1. Compute partials bottom-up (post-order traversal)
        computePartials(tree.getRoot(), siteIndex);

        // 2. Sample root state
        sampleRootState();

        // 3. Sample all descendant states (pre-order traversal)
        sampleDescendantStates(tree.getRoot());

        // 4. Update leaf states in alignment
        updateLeafStates(siteIndex);
    }

    /**
     * Post-order traversal to compute partial likelihoods.
     * For leaves: partial = P(read_counts | genotype)
     * For internal nodes: partial = product over children of sum over child states
     */
    private void computePartials(Node node, int siteIndex) {
        if (node.isLeaf()) {
            computeLeafPartial(node, siteIndex);
        } else {
            // First compute children's partials
            for (int i = 0; i < node.getChildCount(); i++) {
                computePartials(node.getChild(i), siteIndex);
            }
            // Then compute this node's partial
            computeInternalPartial(node);
        }
    }

    /**
     * Compute partial likelihood for a leaf node using read count model.
     * partial[g] = P(read_counts | genotype = g)
     */
    private void computeLeafPartial(Node node, int siteIndex) {
        int nodeNr = node.getNr();
        int taxonIndex = nodeNrToTaxonIndex[nodeNr];
        int[] counts = readCount.getReadCounts(taxonIndex, siteIndex);
        int coverage = coverages[taxonIndex][siteIndex];

        for (int g = 0; g < numStates; g++) {
            double logP = readCountModel.logLiklihoodRC(g, counts, coverage, taxonIndex);
            partials[nodeNr][g] = Math.exp(logP);
        }
    }

    /**
     * Compute partial likelihood for an internal node.
     * partial[g] = product over children of (sum over child states of P(child|g) * partial[child])
     */
    private void computeInternalPartial(Node node) {
        int nodeNr = node.getNr();

        // Initialize partials to 1 (for multiplication)
        Arrays.fill(partials[nodeNr], 1.0);

        // For each child, multiply in the contribution
        for (int childIdx = 0; childIdx < node.getChildCount(); childIdx++) {
            Node child = node.getChild(childIdx);
            int childNr = child.getNr();

            // Get transition matrix for this branch
            getTransitionMatrix(child, transitionMatrix);

            // For each parent state, compute sum over child states
            for (int parentState = 0; parentState < numStates; parentState++) {
                double sum = 0.0;
                for (int childState = 0; childState < numStates; childState++) {
                    // transitionMatrix[parentState * numStates + childState] = P(childState | parentState)
                    double transProb = transitionMatrix[parentState * numStates + childState];
                    sum += transProb * partials[childNr][childState];
                }
                partials[nodeNr][parentState] *= sum;
            }
        }
    }

    /**
     * Get the transition probability matrix for a branch.
     */
    private void getTransitionMatrix(Node node, double[] matrix) {
        double branchRate = (branchRateModel != null) ?
                branchRateModel.getRateForBranch(node) : 1.0;
        double rate = branchRate * siteModel.getRateForCategory(0, node);

        substitutionModel.getTransitionProbabilities(
                node,
                node.getParent().getHeight(),
                node.getHeight(),
                rate,
                matrix
        );
    }

    /**
     * Sample the root state from P(root = g) proportional to pi[g] * partial[root][g]
     */
    private void sampleRootState() {
        Node root = tree.getRoot();
        int rootNr = root.getNr();

        // Get equilibrium frequencies
        rootFrequencies = substitutionModel.getFrequencies();

        // Compute log probabilities
        double maxLogProb = Double.NEGATIVE_INFINITY;
        for (int g = 0; g < numStates; g++) {
            logProbs[g] = Math.log(rootFrequencies[g]) + Math.log(partials[rootNr][g]);
            maxLogProb = Math.max(maxLogProb, logProbs[g]);
        }

        // Sample from normalized distribution
        sampledStates[rootNr] = sampleFromLogProbs(logProbs, maxLogProb);
    }

    /**
     * Pre-order traversal to sample all descendant states.
     * For each node: P(node = g | parent) proportional to P(g | parentState) * partial[node][g]
     */
    private void sampleDescendantStates(Node node) {
        for (int childIdx = 0; childIdx < node.getChildCount(); childIdx++) {
            Node child = node.getChild(childIdx);
            int childNr = child.getNr();
            int parentState = sampledStates[node.getNr()];

            // Get transition matrix
            getTransitionMatrix(child, transitionMatrix);

            // Compute log probabilities for each child state
            double maxLogProb = Double.NEGATIVE_INFINITY;
            for (int g = 0; g < numStates; g++) {
                double transProb = transitionMatrix[parentState * numStates + g];
                logProbs[g] = Math.log(transProb) + Math.log(partials[childNr][g]);
                maxLogProb = Math.max(maxLogProb, logProbs[g]);
            }

            // Sample child state
            sampledStates[childNr] = sampleFromLogProbs(logProbs, maxLogProb);

            // Recurse to children
            if (!child.isLeaf()) {
                sampleDescendantStates(child);
            }
        }
    }

    /**
     * Update leaf states in the mutable alignment.
     */
    private void updateLeafStates(int siteIndex) {
        for (Node leaf : tree.getExternalNodes()) {
            int nodeNr = leaf.getNr();
            int taxonIndex = nodeNrToTaxonIndex[nodeNr];
            leafStates[taxonIndex] = sampledStates[nodeNr];
        }
        alignment.setSiteValuesBySite(siteIndex, leafStates);
    }

    /**
     * Sample a state from log probabilities using log-sum-exp for numerical stability.
     */
    private int sampleFromLogProbs(double[] logProbs, double maxLogProb) {
        // Convert to probabilities using log-sum-exp trick
        double sumExp = 0.0;
        double[] probs = new double[numStates];
        for (int i = 0; i < numStates; i++) {
            probs[i] = Math.exp(logProbs[i] - maxLogProb);
            sumExp += probs[i];
        }

        // Normalize and sample
        double rand = Randomizer.nextDouble() * sumExp;
        double cumulative = 0.0;
        for (int i = 0; i < numStates; i++) {
            cumulative += probs[i];
            if (rand < cumulative) {
                return i;
            }
        }

        // Fallback (should not reach here)
        return numStates - 1;
    }
}
