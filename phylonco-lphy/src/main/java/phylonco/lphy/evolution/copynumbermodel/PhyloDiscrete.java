package phylonco.lphy.evolution.copynumbermodel;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;
import java.util.Map;

public class PhyloDiscrete implements GenerativeDistribution<IntegerCharacterMatrix> {

    private Value<MarkovTraitEvolution<Integer>> model;
    private Value<TimeTree> tree;
    private Value<Integer> length;
    private int[][] data;

    // Parameter names
    public static final String modelParam = "evolutionModel";
    public static final String treeParam = "tree";
    public static final String lengthParam = "length";

    public PhyloDiscrete(
            @ParameterInfo(name = modelParam, description = "The trait evolution model") Value<MarkovTraitEvolution<Integer>> evolutionModel,
            @ParameterInfo(name = treeParam, description = "A tree with taxa") Value<TimeTree> tree,
            @ParameterInfo(name = lengthParam, description = "The number of traits/bins to evolve") Value<Integer> length) {

        super();
        this.model = evolutionModel;
        this.tree = tree;
        this.length = length;
        // Initialize data structure
        int numNodes = tree.value().nodeCount();
        this.data = new int[numNodes][length.value()];
       }

    // Initializes the root state.
    private void initializeRoot(TimeTreeNode rootNode) {
        int rootState = model.value().sampleAncestralTrait();
        for (int i = 0; i < length.value(); i++) {
            data[rootNode.getIndex()][i] = rootState;
        }
    }

    private void evolveBranch(TimeTreeNode parentNode, TimeTreeNode childNode,
                              double startTime, double endTime, IntegerCharacterMatrix matrix) {

        int parentIndex = parentNode.getIndex();
        int childIndex = childNode.getIndex();
        double timeInterval = startTime - endTime; //Branch length

        // For each bin/trait, evolve from parent to child
        for (int i = 0; i < length.value(); i++) {
            int parentState = data[parentIndex][i];
            int childState = model.value().evolveTraitOverTime(parentState, timeInterval);
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
        // Get taxa from the tree
        Taxa taxa = tree.value().getTaxa();

        // Create the output matrix
        IntegerCharacterMatrix matrix = new IntegerCharacterMatrix(taxa, length.value());

        // Simulate trait evolution by traversing the tree
        TimeTreeNode root = tree.value().getRoot();
        traverse(root, matrix);

        // Return the result as a random variable
        return new RandomVariable<>("cnvMatrix", matrix, this);
    }

    @Override
    public Map<String, Value> getParams() {
        return Map.of(
                modelParam, model,
                treeParam, tree,
                lengthParam, length
        );
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(modelParam)) {
            model = value;
        } else if (paramName.equals(treeParam)) {
            tree = value;
        } else if (paramName.equals(lengthParam)) {
            length = value;
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
        return length;
    }
}

