package phylonco.beast.evolution.likelihood;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * Logs a tree-aware sample of the per-cell genotype alignment when genotypes are integrated out by
 * {@link ReadCountTreeLikelihood} (so there is no sampled MutableAlignment to log).
 *
 * <p>For each logged MCMC state it draws a joint genotype sample over the whole tree, exactly as
 * {@code GibbsSiteOperator} (sampleAllSites) did: leaf partials are read-count likelihoods
 * {@code P(reads | g)}, a Felsenstein up-pass computes internal partials per rate category, then a
 * rate category and the root state are sampled and a pre-order down-pass samples every node's
 * genotype. Sampling the ancestral (internal) states is what conditions the leaf genotypes on the
 * tree; only the leaf genotypes are emitted here (one taxon-named column each), the direct analog
 * of the old MutableAlignment log but now conditioned on the tree rather than read counts alone.</p>
 */
@Description("Tree-aware sample of per-cell genotypes (leaf states from a joint ancestral down-pass)")
public class SampledGenotypeLogger extends BEASTObject implements Loggable {

    final public Input<Tree> treeInput = new Input<>(
            "tree", "phylogenetic tree", Input.Validate.REQUIRED);
    final public Input<SiteModel.Base> siteModelInput = new Input<>(
            "siteModel", "site model containing the genotype substitution model", Input.Validate.REQUIRED);
    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>(
            "branchRateModel", "branch rate model for relaxed clock (optional)", Input.Validate.OPTIONAL);
    final public Input<LikelihoodReadCountModel> readCountModelInput = new Input<>(
            "readCountModel", "read-count model supplying P(reads | genotype)", Input.Validate.REQUIRED);
    final public Input<ReadCount> readCountInput = new Input<>(
            "readCount", "read count data", Input.Validate.REQUIRED);

    private Tree tree;
    private SiteModel.Base siteModel;
    private SubstitutionModel substitutionModel;
    private BranchRateModel.Base branchRateModel;
    private LikelihoodReadCountModel readCountModel;
    private ReadCount readCount;
    private DataType datatype;

    private int numNodes;
    private int numStates;
    private int numSites;
    private int numCategories;

    private double[][][] partials;      // [category][node][state]
    private int[] sampledStates;        // [node]
    private double[] transitionMatrix;  // [state*state]
    private double[] logProbs;          // [state]
    private double[] categoryLogProbs;  // [category]
    private double[] rootFrequencies;

    // tree leaf nodeNr -> ReadCount row index
    private int[] nodeNrToRCIndex;
    // ReadCount row order for stable output: taxon name per column
    private String[] rcTaxaNames;
    // accumulates the sampled genotype string for each ReadCount taxon across sites
    private char[][] sampledSeq;        // [rcTaxon][site]

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        siteModel = siteModelInput.get();
        substitutionModel = siteModel.getSubstitutionModel();
        branchRateModel = branchRateModelInput.get();
        readCountModel = readCountModelInput.get();
        readCount = readCountInput.get();
        datatype = readCountModel.getDataType();

        numStates = datatype.getStateCount();
        numSites = readCount.getSiteNumber();
        numNodes = tree.getNodeCount();
        numCategories = siteModel.getCategoryCount();

        rcTaxaNames = readCount.getTaxaNames();

        // map each tree leaf (by taxon name) to its ReadCount row
        nodeNrToRCIndex = new int[numNodes];
        Arrays.fill(nodeNrToRCIndex, -1);
        for (Node leaf : tree.getExternalNodes()) {
            String name = leaf.getID();
            int rcIdx = -1;
            for (int j = 0; j < rcTaxaNames.length; j++) {
                if (name.equals(rcTaxaNames[j].trim())) {
                    rcIdx = j;
                    break;
                }
            }
            if (rcIdx == -1) {
                throw new IllegalArgumentException("tree taxon " + name + " not found in ReadCount data");
            }
            nodeNrToRCIndex[leaf.getNr()] = rcIdx;
        }

        partials = new double[numCategories][numNodes][numStates];
        sampledStates = new int[numNodes];
        transitionMatrix = new double[numStates * numStates];
        logProbs = new double[numStates];
        categoryLogProbs = new double[numCategories];
        sampledSeq = new char[rcTaxaNames.length][numSites];
    }

    @Override
    public void init(PrintStream out) {
        for (String name : rcTaxaNames) {
            out.print("genotype(" + name.trim() + ")\t");
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        // refresh the read-count caches from the current parameter values
        readCountModel.initialize();
        for (int site = 0; site < numSites; site++) {
            sampleSite(site);
        }
        for (int t = 0; t < rcTaxaNames.length; t++) {
            out.print(new String(sampledSeq[t]) + "\t");
        }
    }

    /** Joint sample of all node genotypes at one site (ported from GibbsSiteOperator.sampleSite). */
    private void sampleSite(int site) {
        for (Node leaf : tree.getExternalNodes()) {
            computeLeafPartial(leaf, site);
        }
        for (int cat = 0; cat < numCategories; cat++) {
            computeInternalPartials(tree.getRoot(), cat);
        }
        int category = sampleCategoryAndRootState();
        sampleDescendantStates(tree.getRoot(), category);

        // record the leaf genotypes for this site
        for (Node leaf : tree.getExternalNodes()) {
            int rcIdx = nodeNrToRCIndex[leaf.getNr()];
            sampledSeq[rcIdx][site] = datatype.getCharacter(sampledStates[leaf.getNr()]).charAt(0);
        }
    }

    private void computeLeafPartial(Node node, int site) {
        int nodeNr = node.getNr();
        int rcIdx = nodeNrToRCIndex[nodeNr];
        int[] counts = readCount.getReadCounts(rcIdx, site);
        int coverage = 0;
        for (int c : counts) coverage += c;
        for (int g = 0; g < numStates; g++) {
            double partial = Math.exp(readCountModel.logLiklihoodRC(g, counts, coverage, rcIdx));
            for (int cat = 0; cat < numCategories; cat++) {
                partials[cat][nodeNr][g] = partial;
            }
        }
    }

    private void computeInternalPartials(Node node, int category) {
        if (node.isLeaf()) return;
        for (int i = 0; i < node.getChildCount(); i++) {
            computeInternalPartials(node.getChild(i), category);
        }
        int nodeNr = node.getNr();
        Arrays.fill(partials[category][nodeNr], 1.0);
        for (int childIdx = 0; childIdx < node.getChildCount(); childIdx++) {
            Node child = node.getChild(childIdx);
            int childNr = child.getNr();
            getTransitionMatrix(child, category, transitionMatrix);
            for (int parentState = 0; parentState < numStates; parentState++) {
                double sum = 0.0;
                for (int childState = 0; childState < numStates; childState++) {
                    sum += transitionMatrix[parentState * numStates + childState]
                            * partials[category][childNr][childState];
                }
                partials[category][nodeNr][parentState] *= sum;
            }
        }
    }

    private void getTransitionMatrix(Node node, int category, double[] matrix) {
        double branchRate = (branchRateModel != null) ? branchRateModel.getRateForBranch(node) : 1.0;
        double rate = branchRate * siteModel.getRateForCategory(category, node);
        substitutionModel.getTransitionProbabilities(
                node, node.getParent().getHeight(), node.getHeight(), rate, matrix);
    }

    private int sampleCategoryAndRootState() {
        Node root = tree.getRoot();
        int rootNr = root.getNr();
        rootFrequencies = substitutionModel.getFrequencies();
        double[] proportions = siteModel.getCategoryProportions(root);

        double maxCatLogProb = Double.NEGATIVE_INFINITY;
        for (int cat = 0; cat < numCategories; cat++) {
            double marginal = 0.0;
            for (int g = 0; g < numStates; g++) {
                marginal += rootFrequencies[g] * partials[cat][rootNr][g];
            }
            categoryLogProbs[cat] = Math.log(proportions[cat]) + Math.log(marginal);
            maxCatLogProb = Math.max(maxCatLogProb, categoryLogProbs[cat]);
        }
        int sampledCategory = sampleFromLogProbs(categoryLogProbs, numCategories, maxCatLogProb);

        double maxLogProb = Double.NEGATIVE_INFINITY;
        for (int g = 0; g < numStates; g++) {
            logProbs[g] = Math.log(rootFrequencies[g]) + Math.log(partials[sampledCategory][rootNr][g]);
            maxLogProb = Math.max(maxLogProb, logProbs[g]);
        }
        sampledStates[rootNr] = sampleFromLogProbs(logProbs, numStates, maxLogProb);
        return sampledCategory;
    }

    private void sampleDescendantStates(Node node, int category) {
        for (int childIdx = 0; childIdx < node.getChildCount(); childIdx++) {
            Node child = node.getChild(childIdx);
            int childNr = child.getNr();
            int parentState = sampledStates[node.getNr()];
            getTransitionMatrix(child, category, transitionMatrix);
            double maxLogProb = Double.NEGATIVE_INFINITY;
            for (int g = 0; g < numStates; g++) {
                logProbs[g] = Math.log(transitionMatrix[parentState * numStates + g])
                        + Math.log(partials[category][childNr][g]);
                maxLogProb = Math.max(maxLogProb, logProbs[g]);
            }
            sampledStates[childNr] = sampleFromLogProbs(logProbs, numStates, maxLogProb);
            if (!child.isLeaf()) {
                sampleDescendantStates(child, category);
            }
        }
    }

    private int sampleFromLogProbs(double[] logProbsArr, int size, double maxLogProb) {
        double sumExp = 0.0;
        double[] probs = new double[size];
        for (int i = 0; i < size; i++) {
            probs[i] = Math.exp(logProbsArr[i] - maxLogProb);
            sumExp += probs[i];
        }
        double rand = Randomizer.nextDouble() * sumExp;
        double cumulative = 0.0;
        for (int i = 0; i < size; i++) {
            cumulative += probs[i];
            if (rand < cumulative) return i;
        }
        return size - 1;
    }

    @Override
    public void close(PrintStream out) {
    }
}
