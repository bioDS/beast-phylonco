package phylonco.lphy.evolution.copynumbermodel;

import lphy.base.evolution.alignment.TaxaCharacterMatrix;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class NestedBD implements GenerativeDistribution<TaxaCharacterMatrix<Integer>> {

    private final Value<Integer> nBins;
    private final Value<Integer> nCells;
    private final Value<Double> lambda;
    private final Value<Double> mu;
    private final Value<TimeTree> tree;
    private final int[][] data;

    public static final String birthRateString = "lambda";
    public static final String deathRateString = "mu";
    public static final String nBinsString = "nBins";
    public static final String nCellsString = "nCells";
    public static final String treeString = "tree";

    public NestedBD(
            @ParameterInfo(name = birthRateString, description = "The birth rate (generation) of copy number variants") Value<Double> lambda,
            @ParameterInfo(name = deathRateString, description = "The death rate (loss) of copy number variants") Value<Double> mu,
            @ParameterInfo(name = nBinsString, description = "Number of genomic bins") Value<Integer> nBins,
            @ParameterInfo(name = nCellsString, description = "Number of cells") Value<Integer> nCells,
            @ParameterInfo(name = treeString, description = "Number ") Value<TimeTree> tree

    ) {
        super();
        this.lambda = lambda;
        this.mu = mu;
        this.nBins = nBins;
        this.nCells = nCells;
        this.tree = tree;

        int numNodes = tree.value().nodeCount(); // total internal nodes and external nodes
        this.data = new int[nBins.value()][numNodes]; // set up bins and node data matrix

        // the root has 2 copies
        for (int i = 0; i < nBins.value(); i++) {
            this.data[i][tree.value().getRoot().getIndex()] = 2; // 2 copies
        }
    }


    private void traverse(TimeTreeNode node) {

        double tEnd = node.getAge();
        double tStart = node.getParent().getAge();
        int parentIndex = node.getParent().getIndex(); //Not sure why we need to getParent() twice?
        int currentIndex = node.getIndex();  // Get current node's index

        // Special case for root node
        if (node.isRoot()) {
            // Root already has 2 copies initialised in the constructor
            // Just traverse its children
            traverse(node.getLeft());
            traverse(node.getRight());
            return;
        }

        if (node.isLeaf()) {
            // simulate down branch parent to leaf
            for (int i = 0; i < nBins.value(); i++) {
                // generate copy number for child from parent
                simulateCopiesOnBranch(tStart, tEnd, this.data[i][parentIndex], currentIndex);
            }
        } else {
            // simulate for one branch
            // for each bin
            for (int i = 0; i < nBins.value(); i++) {
                // generate copy number for child from parent
                simulateCopiesOnBranch(tStart, tEnd, this.data[i][parentIndex], currentIndex);
            }

            // repeats for left and right children
            traverse(node.getLeft()); // left child
            traverse(node.getRight()); // right child
        }
    }

    private void simulateCopiesOnBranch(double startTime, double endTime, int startCopies, int nodeIndex) {
        // For each bin
        for (int binIdx = 0; binIdx < nBins.value(); binIdx++) {
            double tCurrent = startTime;
            int m = startCopies; // Current copy number
            // Simulate events until we reach the end time or copies go extinct
            while (tCurrent < endTime && m > 0) {
                // Calculate event rates
                double birthRate = lambda.value() * m;
                double deathRate = mu.value() * m;
                double totalRate = birthRate + deathRate;

                // If no events are possible, break
                if (totalRate == 0) {
                    break;
                }

                // Sample time to next event from exponential distribution
                double tNext = -Math.log(Math.random()) / totalRate;

                // Check if we've reached branch end
                if (tCurrent + tNext > endTime) {
                    break;
                }
                // Update current time
                tCurrent += tNext;
                // Determine event type
                double u = Math.random();
                if (u < birthRate / totalRate) {
                    // Birth event (copy gain)
                    m += 1;
                } else {
                    // Death event (copy loss)
                    m -= 1;
                }
            }
            // Update copy number for this bin
            this.data[binIdx][nodeIndex] = m;
        }

    }


    @Override
    public RandomVariable<TaxaCharacterMatrix<Integer>> sample() {
        // sample your data bins x cells matrix
        // Get parameter values
        int numBins = nBins.value();
        int numCells = nCells.value();
        double lambdaVal = lambda.value();
        double muVal = mu.value();


        // Need a random number generator
        double randomDouble = Math.random();

        //Tree topology
        TimeTree treeValue = tree.value();
        TimeTreeNode root = treeValue.getRoot();
        // traverse tree by breadth first search

        // example doing first bin
        int binIndex = 0;
        traverse(root);

        return null;
    }

    @Override
    public Map<String, Value> getParams() {
        return Map.of();
    }
}
