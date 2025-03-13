package phylonco.lphy.evolution.copynumbermodel;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class NestedBD implements GenerativeDistribution<CopyNumberProfiles> {

    private Value<Integer> nBins;
    private Value<Integer> nCells;
    private Value<Double> lambda;
    private Value<Double> mu;
    private Value<TimeTree> tree;
    private int[][] data;

    private CopyNumberProfiles copyProfiles;

    public static final String birthRateString = "lambda";
    public static final String deathRateString = "mu";
    public static final String nBinsString = "nBins";
    public static final String nCellsString = "nCells";
    public static final String treeString = "tree";

    // for testing
    public NestedBD() {
    }

    public NestedBD(
            @ParameterInfo(name = birthRateString, description = "The birth rate (generation) of copy number variants") Value<Double> lambda,
            @ParameterInfo(name = deathRateString, description = "The death rate (loss) of copy number variants") Value<Double> mu,
            @ParameterInfo(name = nBinsString, description = "Number of genomic bins") Value<Integer> nBins,
            @ParameterInfo(name = nCellsString, description = "Number of cells") Value<Integer> nCells,
            @ParameterInfo(name = treeString, description = "Tree with taxa") Value<TimeTree> tree

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

        // Special case for root node
        if (node.isRoot()) {
            // Root already has 2 copies initialised in the constructor
            // Just traverse its children
            traverse(node.getLeft());
            traverse(node.getRight());
            return;
        }

        // Get the parent and child's age
        double tEnd = node.getAge();
        double tStart = node.getParent().getAge();
        int parentIndex = node.getParent().getIndex();

        if (node.isLeaf()) {
            // simulate down branch parent to leaf
            for (int i = 0; i < nBins.value(); i++) {
                // generate copy number for child from parent
                simulateCopiesOnBranch(tStart, tEnd, this.data[i][parentIndex], node);
                // put data into CopyNumberProfile
            }

        } else {
            // simulate for one branch
            // for each bin
            for (int i = 0; i < nBins.value(); i++) {
                // generate copy number for child from parent
                simulateCopiesOnBranch(tStart, tEnd, this.data[i][parentIndex], node);
            }

            // repeats for left and right children
            traverse(node.getLeft()); // left child
            traverse(node.getRight()); // right child
        }
    }

    // single bin
    public int simulateCopiesOnBranchBin(double startTime, double endTime, int startCopies) {
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
            RandomGenerator randomGen = RandomUtils.getRandom();
            double tNext = -Math.log(randomGen.nextDouble()) / totalRate;

            // Check if we've reached branch end
            if (tCurrent + tNext > endTime) {
                break;
            }
            // Update current time
            tCurrent += tNext;
            // Determine event type
            double u = randomGen.nextDouble();
            if (u < birthRate / totalRate) {
                // Birth event (copy gain)
                m += 1;
            } else {
                // Death event (copy loss)
                m -= 1;
            }
        }
        return m;

    }

    // all bins
    private void simulateCopiesOnBranch(double startTime, double endTime, int startCopies, TimeTreeNode node) {
        // For each bin
        for (int binIdx = 0; binIdx < nBins.value(); binIdx++) {
            // Update copy number for this bin
            int m = simulateCopiesOnBranchBin(startTime, endTime, startCopies);
            this.data[binIdx][node.getIndex()] = m;
            // if node is a leaf put data into copy number profile
            if (node.isLeaf()) {
                this.copyProfiles.setState(node.getIndex(), binIdx, m);
            }
        }

    }


    @Override
    public RandomVariable<CopyNumberProfiles> sample() {
        // get Taxa from the tree
        Taxa taxa = tree.value().getTaxa();
        this.copyProfiles = new CopyNumberProfiles(taxa, this.nBins.value());

        //Tree topology
        TimeTree treeValue = tree.value();
        TimeTreeNode root = treeValue.getRoot();
        // traverse tree by breadth first search

        // example doing first bin
        int binIndex = 0;
        traverse(root);

        // return that copy number profile

//        return new RandomVariable(copyProfiles;

        return new RandomVariable<CopyNumberProfiles>(null, copyProfiles, this);
    }

    @Override
    public Map<String, Value> getParams() {
        return Map.of(
                birthRateString, lambda,
                deathRateString, mu,
                nBinsString, nBins,
                nCellsString, nCells,
                treeString, tree
        );
    }


    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(birthRateString)) {
            lambda = value;
        } else if (paramName.equals(deathRateString)) {
            mu = value;
        } else if (paramName.equals(nBinsString)) {
            nBins = value;
        } else if (paramName.equals(nCellsString)) {
            nCells = value;
        } else if (paramName.equals(treeString)) {
            tree = value;
        } else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }

    // getParameter() for each parameter that returns Value<Type>
    public Value<Double> getLambda() {
        return lambda;
    }

    public Value<Double> getMu() {
        return mu;
    }

    public Value<Integer> getNBins() {
        return nBins;
    }

    public Value<Integer> getNCells() {
        return nCells;
    }

    public Value<TimeTree> getTree() {
        return tree;
    }
}
package phylonco.lphy.evolution.copynumbermodel;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class NestedBD implements GenerativeDistribution<CopyNumberProfiles> {

    private Value<Integer> nBins;
    private Value<Integer> nCells;
    private Value<Double> lambda;
    private Value<Double> mu;
    private Value<TimeTree> tree;
    private int[][] data;

    private CopyNumberProfiles copyProfiles;

    public static final String birthRateString = "lambda";
    public static final String deathRateString = "mu";
    public static final String nBinsString = "nBins";
    public static final String nCellsString = "nCells";
    public static final String treeString = "tree";

    // for testing
    public NestedBD() {
    }

    public NestedBD(
            @ParameterInfo(name = birthRateString, description = "The birth rate (generation) of copy number variants") Value<Double> lambda,
            @ParameterInfo(name = deathRateString, description = "The death rate (loss) of copy number variants") Value<Double> mu,
            @ParameterInfo(name = nBinsString, description = "Number of genomic bins") Value<Integer> nBins,
            @ParameterInfo(name = nCellsString, description = "Number of cells") Value<Integer> nCells,
            @ParameterInfo(name = treeString, description = "Tree with taxa") Value<TimeTree> tree

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

        // Special case for root node
        if (node.isRoot()) {
            // Root already has 2 copies initialised in the constructor
            // Just traverse its children
            traverse(node.getLeft());
            traverse(node.getRight());
            return;
        }

        // Get the parent and child's age
        double tEnd = node.getAge();
        double tStart = node.getParent().getAge();
        int parentIndex = node.getParent().getIndex();

        if (node.isLeaf()) {
            // simulate down branch parent to leaf
            for (int i = 0; i < nBins.value(); i++) {
                // generate copy number for child from parent
                simulateCopiesOnBranch(tStart, tEnd, this.data[i][parentIndex], node);
                // put data into CopyNumberProfile
            }

        } else {
            // simulate for one branch
            // for each bin
            for (int i = 0; i < nBins.value(); i++) {
                // generate copy number for child from parent
                simulateCopiesOnBranch(tStart, tEnd, this.data[i][parentIndex], node);
            }

            // repeats for left and right children
            traverse(node.getLeft()); // left child
            traverse(node.getRight()); // right child
        }
    }

    // single bin
    public int simulateCopiesOnBranchBin(double startTime, double endTime, int startCopies) {
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
            RandomGenerator randomGen = RandomUtils.getRandom();
            double tNext = -Math.log(randomGen.nextDouble()) / totalRate;

            // Check if we've reached branch end
            if (tCurrent + tNext > endTime) {
                break;
            }
            // Update current time
            tCurrent += tNext;
            // Determine event type
            double u = randomGen.nextDouble();
            if (u < birthRate / totalRate) {
                // Birth event (copy gain)
                m += 1;
            } else {
                // Death event (copy loss)
                m -= 1;
            }
        }
        return m;

    }

    // all bins
    private void simulateCopiesOnBranch(double startTime, double endTime, int startCopies, TimeTreeNode node) {
        // For each bin
        for (int binIdx = 0; binIdx < nBins.value(); binIdx++) {
            // Update copy number for this bin
            int m = simulateCopiesOnBranchBin(startTime, endTime, startCopies);
            this.data[binIdx][node.getIndex()] = m;
            // if node is a leaf put data into copy number profile
            if (node.isLeaf()) {
                this.copyProfiles.setState(node.getIndex(), binIdx, m);
            }
        }

    }


    @Override
    public RandomVariable<CopyNumberProfiles> sample() {
        // get Taxa from the tree
        Taxa taxa = tree.value().getTaxa();
        this.copyProfiles = new CopyNumberProfiles(taxa, this.nBins.value());

        //Tree topology
        TimeTree treeValue = tree.value();
        TimeTreeNode root = treeValue.getRoot();
        // traverse tree by breadth first search

        // example doing first bin
        int binIndex = 0;
        traverse(root);

        // return that copy number profile

//        return new RandomVariable(copyProfiles;

        return new RandomVariable<CopyNumberProfiles>(null, copyProfiles, this);
    }

    @Override
    public Map<String, Value> getParams() {
        return Map.of(
                birthRateString, lambda,
                deathRateString, mu,
                nBinsString, nBins,
                nCellsString, nCells,
                treeString, tree
        );
    }


    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(birthRateString)) {
            lambda = value;
        } else if (paramName.equals(deathRateString)) {
            mu = value;
        } else if (paramName.equals(nBinsString)) {
            nBins = value;
        } else if (paramName.equals(nCellsString)) {
            nCells = value;
        } else if (paramName.equals(treeString)) {
            tree = value;
        } else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }

    // getParameter() for each parameter that returns Value<Type>
    public Value<Double> getLambda() {
        return lambda;
    }

    public Value<Double> getMu() {
        return mu;
    }

    public Value<Integer> getNBins() {
        return nBins;
    }

    public Value<Integer> getNCells() {
        return nCells;
    }

    public Value<TimeTree> getTree() {
        return tree;
    }
}
