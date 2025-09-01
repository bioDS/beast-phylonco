package phylonco.lphy.evolution.copynumbermodel;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

public class PhyloDiscrete implements GenerativeDistribution<IntegerCharacterMatrix> {

    private Value<MarkovTraitEvolution<Integer>> model;
    private Value<TimeTree> tree;
    private Value<Integer> L;
    private int[][] data;
    private Value<Number> clockRate;
    private Value<Double[]> branchRates;
    private Value<IntegerCharacterMatrix> realData;

    // Parameter names
    public static final String modelParam = "evolutionModel";
    public static final String treeParam = "tree";
    public static final String lengthParam = "L";
    public static final String muParam = "mu";
    public static final String branchRatesParam = "branchRates";
    public static final String realDataParam = "realData";

    public PhyloDiscrete(
            @ParameterInfo(name = modelParam, description = "The trait evolution model", optional = true) Value<MarkovTraitEvolution<Integer>> evolutionModel,
            @ParameterInfo(name = treeParam, description = "A tree with taxa") Value<TimeTree> tree,
            @ParameterInfo(name = lengthParam, description = "The number of traits/bins to evolve") Value<Integer> L,
            @ParameterInfo(name = muParam, narrativeName = "molecular clock rate", description = "the clock rate. Default value is 1.0.", optional = true)
            Value<Number> mu,
            @ParameterInfo(name = branchRatesParam, description = "a rate for each branch in the tree. Branch rates are assumed to be 1.0 otherwise.", optional = true)
            Value<Double[]> branchRates,
            @ParameterInfo(name = realDataParam, description = "Real CNV data matrix to use instead of simulation", optional = true) Value<IntegerCharacterMatrix> realData){

        super();
        // Validation: either realData OR (evolutionModel + L) must be provided
        if (realData == null && (evolutionModel == null || L == null)) {
            throw new IllegalArgumentException("Either realData must be provided, or both evolutionModel and L must be provided for simulation");
        }
        this.model = evolutionModel;
        this.tree = tree;
        this.L = L;
        this.clockRate = mu;
        this.realData = realData;

        Double[] treeBranchRates = tree.value().getBranchRates();

        if (treeBranchRates != null && treeBranchRates.length > 0) {
            if (branchRates != null) {
                LoggerUtils.log.warning("PhyloDiscrete has branchRates from input parameter and tree has branch rates, " +
                        "default to using input parameter branchRates.");
            } else {
                this.branchRates = new Value<>("branchRates", treeBranchRates);
            }
        }
        this.branchRates = branchRates;

        // Initialize data structure for simulation only
        if (realData == null) {
            int numNodes = tree.value().nodeCount();
            this.data = new int[numNodes][L.value()];
        }
    }

    // Initializes the root state.
    private void initializeRoot(TimeTreeNode rootNode) {
        int rootState = model.value().sampleAncestralTrait();
        for (int i = 0; i < L.value(); i++) {
            data[rootNode.getIndex()][i] = rootState;
        }
    }

    private void evolveBranch(TimeTreeNode parentNode, TimeTreeNode childNode,
                              double startTime, double endTime, IntegerCharacterMatrix matrix) {

        int parentIndex = parentNode.getIndex();
        int childIndex = childNode.getIndex();
        double timeInterval = startTime - endTime; //Branch length
        double mu = (this.clockRate == null) ? 1.0 : ValueUtils.doubleValue(clockRate); //clock rate default = 1.0
        double evolutionaryBranchLength = timeInterval * mu;

        // Apply branch-specific rate if available
        if (branchRates != null) {
            double branchRateMultiplier = branchRates.value()[childIndex];
            evolutionaryBranchLength *= branchRateMultiplier;
        }

        // For each bin/trait, evolve from parent to child
        for (int i = 0; i < L.value(); i++) {
            int parentState = data[parentIndex][i];
            int childState = model.value().evolveTraitOverTime(parentState, evolutionaryBranchLength);
//            System.out.println("Bin " + i + ": Evolving from parent state " + parentState); //debug

            // Store result
            data[childIndex][i] = childState;

            // If child is a leaf, update the output matrix
            if (childNode.isLeaf()) {
                matrix.setState(childNode.getIndex(), i, childState);
            }
        }
    }

    // Traverses the tree to simulate trait evolution.
    private void traverse(TimeTreeNode node, IntegerCharacterMatrix matrix) {
        // Special case for the root
        if (node.isRoot()) {
            initializeRoot(node);

            // Traverse children
            traverse(node.getLeft(), matrix);
            traverse(node.getRight(), matrix);
            return;
        }

        // Get branch times
        double endTime = node.getAge();
        double startTime = node.getParent().getAge();

        // Evolve traits along this branch
        evolveBranch(node.getParent(), node, startTime, endTime, matrix);

        // Continue traversal if not a leaf
        if (!node.isLeaf()) {
            traverse(node.getLeft(), matrix);
            traverse(node.getRight(), matrix);
        }
    }

    @Override
    public RandomVariable<IntegerCharacterMatrix> sample() {
        // If real data is provided, use it directly
        if (realData != null) {
            return new RandomVariable<>("cnvMatrix", realData.value(), this);
        }

        // Otherwise, Get simulated taxa from the tree
        Taxa taxa = tree.value().getTaxa();

        // Create the output matrix
        IntegerCharacterMatrix matrix = new IntegerCharacterMatrix(taxa, L.value());

        // Simulate trait evolution by traversing the tree
        TimeTreeNode root = tree.value().getRoot();
        traverse(root, matrix);

        // Return the result as a random variable
        return new RandomVariable<>("cnvMatrix", matrix, this);
    }

    @Override
    public Map<String, Value> getParams() {
        SortedMap<String, Value> map = new TreeMap<>();
        map.put(modelParam, model);
        map.put(treeParam, tree);
        map.put(lengthParam, L);
        if (clockRate != null) map.put(muParam, clockRate);
        if (branchRates != null) map.put(branchRatesParam, branchRates);
        if (realData != null) map.put(realDataParam, realData);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(modelParam)) {
            model = value;
        } else if (paramName.equals(treeParam)) {
            tree = value;
        } else if (paramName.equals(lengthParam)) {
            L = value;
        } else if (paramName.equals(muParam)) {
            clockRate = value;
        } else if (paramName.equals(branchRatesParam)) {
            branchRates = value;
        } else if (paramName.equals(realDataParam)) {
            realData = value;
        } else {
            throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }
    }

    // Getters
    public Value<MarkovTraitEvolution<Integer>> getModel() {
        return model;
    }

    public Value<TimeTree> getTree() {
        return tree;
    }

    public Value<Integer> getLength() {
        return L;
    }

    public Value<Number> getClockRate() {
        return clockRate;
    }

    public Value<Double[]> getBranchRates() {
        return branchRates;
    }
}

