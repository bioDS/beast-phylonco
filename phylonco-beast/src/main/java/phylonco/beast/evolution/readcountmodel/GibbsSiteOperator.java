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
import phylonco.beast.evolution.datatype.ReadCountMatrix;

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
 * Integrates over site rate categories: for each site, partials are computed
 * under every rate category, a category is sampled, and node states are then
 * sampled using that category's transition matrices.
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

    public Input<ReadCountMatrix> readCountInput = new Input<>(
            "readCount",
            "read count data",
            Input.Validate.REQUIRED);

    public Input<Boolean> sampleAllSitesInput = new Input<>(
            "sampleAllSites",
            "if true, sample all sites in one proposal; if false (default), sample sites according to numSitesToResample",
            false);

    public Input<Integer> numSitesToResampleInput = new Input<>(
            "numSitesToResample",
            "number of sites to resample per proposal when sampleAllSites=false. " +
                    "Default 1 (original behaviour). Sites are sampled without replacement.",
            1);

    public Input<Boolean> useEntropyWeightingInput = new Input<>(
            "useEntropyWeighting",
            "weight site-selection probability by per-site pooled-nucleotide Shannon entropy " +
                    "(sum across leaves). Detailed balance preserved (always-accepted Gibbs; " +
                    "site selection probability is a function of data only). Default false.",
            false);

    // Cached references
    private MutableAlignment alignment;
    private Tree tree;
    private SubstitutionModel substitutionModel;
    private SiteModel.Base siteModel;
    private BranchRateModel.Base branchRateModel;
    private LikelihoodReadCountModel readCountModel;
    private ReadCountMatrix readCount;
    private boolean sampleAllSites;
    private int numSitesToResample;
    private boolean useEntropyWeighting;
    // Per-site entropy weight (pooled-nucleotide Shannon entropy), filled lazily on first proposal.
    private double[] siteEntropy;

    // Dimensions
    private int numNodes;
    private int numStates;
    private int numSites;
    private int numTaxa;
    private int numCategories;

    // Pre-allocated arrays for efficiency
    private double[][][] partials;      // [numCategories][numNodes][numStates] - partial likelihoods
    private int[] sampledStates;        // [numNodes] - sampled state for each node
    private int[] leafStates;           // [numTaxa] - leaf states for alignment update
    private double[] transitionMatrix;  // [numStates * numStates] - reusable matrix
    private double[] logProbs;          // [numStates] - reusable probability array
    private double[] categoryLogProbs;  // [numCategories] - for sampling rate category
    private int[][] coverages;          // [numTaxa][numSites] - precomputed coverages
    private double[] rootFrequencies;   // [numStates] - equilibrium frequencies

    // Mapping from node number to taxon index
    private int[] nodeNrToTaxonIndex;
    // Mapping from alignment taxon index to ReadCount taxon index
    private int[] alignToRCIndex;

    // Cached cumulative weights for weighted site sampling
    private int[] cumulativeWeights;
    private int totalWeight;

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
        numSitesToResample = numSitesToResampleInput.get();
        useEntropyWeighting = useEntropyWeightingInput.get();
        if (numSitesToResample < 1) {
            throw new IllegalArgumentException("numSitesToResample must be >= 1 (use sampleAllSites=true for all)");
        }

        numStates = alignment.getDataType().getStateCount();
        numSites = alignment.getSiteCount();
        numTaxa = alignment.getTaxonCount();
        numNodes = tree.getNodeCount();
        numCategories = siteModel.getCategoryCount();

        // Build mapping from alignment taxon index to ReadCount taxon index
        String[] rcTaxaNames = readCount.getTaxonNames();
        alignToRCIndex = new int[numTaxa];
        for (int i = 0; i < numTaxa; i++) {
            String alignTaxonName = alignment.getTaxaNames().get(i);
            alignToRCIndex[i] = i; // default: identity mapping
            for (int j = 0; j < rcTaxaNames.length; j++) {
                if (alignTaxonName.equals(rcTaxaNames[j].trim())) {
                    alignToRCIndex[i] = j;
                    break;
                }
            }
        }

        // Pre-allocate arrays
        partials = new double[numCategories][numNodes][numStates];
        sampledStates = new int[numNodes];
        leafStates = new int[numTaxa];
        transitionMatrix = new double[numStates * numStates];
        logProbs = new double[numStates];
        categoryLogProbs = new double[numCategories];

        // Precompute coverages for all taxa and sites (indexed by alignment taxon order)
        coverages = new int[numTaxa][numSites];
        for (int i = 0; i < numTaxa; i++) {
            for (int j = 0; j < numSites; j++) {
                int[] counts = readCount.getReadCounts(alignToRCIndex[i], j);
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

        // Build cumulative weights for weighted site sampling
        cumulativeWeights = new int[numSites];
        updateCumulativeWeights();
    }

    /**
     * Rebuild cumulative weight array from current pattern weights.
     * Call this after pattern weights change (e.g. after DataTemperedMCMC adds sites).
     */
    public void updateCumulativeWeights() {
        int[] weights = alignment.getWeights();
        totalWeight = 0;
        for (int i = 0; i < numSites; i++) {
            totalWeight += weights[i];
            cumulativeWeights[i] = totalWeight;
        }
    }

    @Override
    public double proposal() {
        int[] weights = alignment.getWeights();
        if (sampleAllSites) {
            // Sample all active sites (skip sites with weight 0)
            for (int siteIndex = 0; siteIndex < numSites; siteIndex++) {
                if (weights[siteIndex] > 0) {
                    sampleSite(siteIndex);
                }
            }
        } else if (numSitesToResample == 1 && !useEntropyWeighting) {
            // Original fast path: sample one site weighted by patternWeight via cumulative search.
            if (totalWeight > 0) {
                int r = Randomizer.nextInt(totalWeight);
                int siteIndex = Arrays.binarySearch(cumulativeWeights, 0, numSites, r + 1);
                if (siteIndex < 0) {
                    siteIndex = -siteIndex - 1;
                }
                sampleSite(siteIndex);
            }
        } else {
            // K-site sampling without replacement, optionally weighted by entropy.
            if (useEntropyWeighting && siteEntropy == null) {
                computeSiteEntropies();
            }
            int[] picked = pickSites(weights, numSitesToResample);
            if (picked != null) {
                for (int s : picked) sampleSite(s);
            }
        }
        // Gibbs operator: always accept
        return Double.POSITIVE_INFINITY;
    }

    /**
     * Pick K distinct sites without replacement, weighted by patternWeight × siteEntropy
     * (if useEntropyWeighting) or patternWeight alone.
     */
    private int[] pickSites(int[] weights, int K) {
        double[] effective = new double[numSites];
        int eligible = 0;
        double totalEffective = 0.0;
        for (int s = 0; s < numSites; s++) {
            if (weights[s] <= 0) continue;
            double w = useEntropyWeighting ? weights[s] * siteEntropy[s] : weights[s];
            if (w > 0.0) {
                effective[s] = w;
                totalEffective += w;
                eligible++;
            }
        }
        if (eligible < K || totalEffective <= 0.0) return null;
        int[] picked = new int[K];
        for (int k = 0; k < K; k++) {
            double r = Randomizer.nextDouble() * totalEffective;
            double cum = 0.0;
            int chosen = -1;
            for (int s = 0; s < numSites; s++) {
                if (effective[s] > 0.0) {
                    cum += effective[s];
                    if (r < cum) { chosen = s; break; }
                }
            }
            if (chosen < 0) {
                for (int s = numSites - 1; s >= 0; s--) {
                    if (effective[s] > 0.0) { chosen = s; break; }
                }
            }
            picked[k] = chosen;
            totalEffective -= effective[chosen];
            effective[chosen] = 0.0;
        }
        return picked;
    }

    /**
     * Compute per-site pooled-nucleotide Shannon entropy across leaves.
     * Same as ExchangeSiteGibbsOperator's siteEntropy. Function of read-count data only,
     * so detailed balance is preserved.
     */
    private void computeSiteEntropies() {
        siteEntropy = new double[numSites];
        double floor = 0.0;
        for (int s = 0; s < numSites; s++) {
            int[] pooled = null;
            int total = 0;
            for (int t = 0; t < numTaxa; t++) {
                int[] counts = readCount.getReadCounts(alignToRCIndex[t], s);
                if (pooled == null) {
                    pooled = new int[counts.length];
                    floor = 0.05 * Math.log(counts.length);  // 5% of max entropy
                }
                for (int n = 0; n < counts.length; n++) {
                    pooled[n] += counts[n];
                    total += counts[n];
                }
            }
            double h = 0.0;
            if (total > 0 && pooled != null) {
                for (int n = 0; n < pooled.length; n++) {
                    if (pooled[n] > 0) {
                        double p = (double) pooled[n] / total;
                        h -= p * Math.log(p);
                    }
                }
            }
            siteEntropy[s] = h + floor;
        }
    }

    /**
     * Sample all leaf genotypes at a single site.
     */
    private void sampleSite(int siteIndex) {
        // 1. Compute leaf partials (independent of rate category)
        for (Node leaf : tree.getExternalNodes()) {
            computeLeafPartial(leaf, siteIndex);
        }

        // 2. For each rate category, compute internal node partials
        for (int cat = 0; cat < numCategories; cat++) {
            computeInternalPartials(tree.getRoot(), cat);
        }

        // 3. Sample rate category and root state
        int category = sampleCategoryAndRootState();

        // 4. Sample all descendant states (pre-order traversal)
        sampleDescendantStates(tree.getRoot(), category);

        // 5. Update leaf states in alignment
        updateLeafStates(siteIndex);
    }

    /**
     * Compute partial likelihood for a leaf node using read count model.
     * Leaf partials are independent of rate category, so the same values
     * are stored across all categories.
     * partial[g] = P(read_counts | genotype = g)
     */
    private void computeLeafPartial(Node node, int siteIndex) {
        int nodeNr = node.getNr();
        int taxonIndex = nodeNrToTaxonIndex[nodeNr];
        int rcIndex = alignToRCIndex[taxonIndex];
        int[] counts = readCount.getReadCounts(rcIndex, siteIndex);
        int coverage = coverages[taxonIndex][siteIndex];

        for (int g = 0; g < numStates; g++) {
            double logP = readCountModel.logLiklihoodRC(g, counts, coverage, rcIndex);
            double partial = Math.exp(logP);
            for (int cat = 0; cat < numCategories; cat++) {
                partials[cat][nodeNr][g] = partial;
            }
        }
    }

    /**
     * Recursively compute internal node partials for a specific rate category.
     * Leaf partials must already be set before calling this.
     */
    private void computeInternalPartials(Node node, int category) {
        if (node.isLeaf()) return;

        for (int i = 0; i < node.getChildCount(); i++) {
            computeInternalPartials(node.getChild(i), category);
        }

        int nodeNr = node.getNr();

        // Initialize partials to 1 (for multiplication)
        Arrays.fill(partials[category][nodeNr], 1.0);

        // For each child, multiply in the contribution
        for (int childIdx = 0; childIdx < node.getChildCount(); childIdx++) {
            Node child = node.getChild(childIdx);
            int childNr = child.getNr();

            // Get transition matrix for this branch under this rate category
            getTransitionMatrix(child, category, transitionMatrix);

            // For each parent state, compute sum over child states
            for (int parentState = 0; parentState < numStates; parentState++) {
                double sum = 0.0;
                for (int childState = 0; childState < numStates; childState++) {
                    double transProb = transitionMatrix[parentState * numStates + childState];
                    sum += transProb * partials[category][childNr][childState];
                }
                partials[category][nodeNr][parentState] *= sum;
            }
        }
    }

    /**
     * Get the transition probability matrix for a branch under a specific rate category.
     */
    private void getTransitionMatrix(Node node, int category, double[] matrix) {
        double branchRate = (branchRateModel != null) ?
                branchRateModel.getRateForBranch(node) : 1.0;
        double rate = branchRate * siteModel.getRateForCategory(category, node);

        substitutionModel.getTransitionProbabilities(
                node,
                node.getParent().getHeight(),
                node.getHeight(),
                rate,
                matrix
        );
    }

    /**
     * Sample the rate category and root state jointly.
     * P(cat, root=g) proportional to proportion[cat] * pi[g] * partials[cat][root][g]
     * First marginalizes over root states to sample the category,
     * then samples the root state within the selected category.
     * Returns the sampled category index.
     */
    private int sampleCategoryAndRootState() {
        Node root = tree.getRoot();
        int rootNr = root.getNr();

        rootFrequencies = substitutionModel.getFrequencies();
        double[] proportions = siteModel.getCategoryProportions(root);

        // Compute marginal log-likelihood for each category
        double maxCatLogProb = Double.NEGATIVE_INFINITY;
        for (int cat = 0; cat < numCategories; cat++) {
            double marginal = 0.0;
            for (int g = 0; g < numStates; g++) {
                marginal += rootFrequencies[g] * partials[cat][rootNr][g];
            }
            categoryLogProbs[cat] = Math.log(proportions[cat]) + Math.log(marginal);
            maxCatLogProb = Math.max(maxCatLogProb, categoryLogProbs[cat]);
        }

        // Sample rate category
        int sampledCategory = sampleFromLogProbs(categoryLogProbs, numCategories, maxCatLogProb);

        // Sample root state within the selected category
        double maxLogProb = Double.NEGATIVE_INFINITY;
        for (int g = 0; g < numStates; g++) {
            logProbs[g] = Math.log(rootFrequencies[g]) + Math.log(partials[sampledCategory][rootNr][g]);
            maxLogProb = Math.max(maxLogProb, logProbs[g]);
        }
        sampledStates[rootNr] = sampleFromLogProbs(logProbs, numStates, maxLogProb);

        return sampledCategory;
    }

    /**
     * Pre-order traversal to sample all descendant states using the selected rate category.
     * For each node: P(node = g | parent) proportional to P(g | parentState) * partial[node][g]
     */
    private void sampleDescendantStates(Node node, int category) {
        for (int childIdx = 0; childIdx < node.getChildCount(); childIdx++) {
            Node child = node.getChild(childIdx);
            int childNr = child.getNr();
            int parentState = sampledStates[node.getNr()];

            // Get transition matrix
            getTransitionMatrix(child, category, transitionMatrix);

            // Compute log probabilities for each child state
            double maxLogProb = Double.NEGATIVE_INFINITY;
            for (int g = 0; g < numStates; g++) {
                double transProb = transitionMatrix[parentState * numStates + g];
                logProbs[g] = Math.log(transProb) + Math.log(partials[category][childNr][g]);
                maxLogProb = Math.max(maxLogProb, logProbs[g]);
            }

            // Sample child state
            sampledStates[childNr] = sampleFromLogProbs(logProbs, numStates, maxLogProb);

            // Recurse to children
            if (!child.isLeaf()) {
                sampleDescendantStates(child, category);
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
    private int sampleFromLogProbs(double[] logProbs, int size, double maxLogProb) {
        // Convert to probabilities using log-sum-exp trick
        double sumExp = 0.0;
        double[] probs = new double[size];
        for (int i = 0; i < size; i++) {
            probs[i] = Math.exp(logProbs[i] - maxLogProb);
            sumExp += probs[i];
        }

        // Normalize and sample
        double rand = Randomizer.nextDouble() * sumExp;
        double cumulative = 0.0;
        for (int i = 0; i < size; i++) {
            cumulative += probs[i];
            if (rand < cumulative) {
                return i;
            }
        }

        // Fallback (should not reach here)
        return size - 1;
    }
}