package phylonco.beast.evolution.tools;

import beast.base.core.BEASTInterface;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.MCMC;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.parser.NexusParser;
import beast.base.parser.XMLParser;
import beast.base.spec.inference.parameter.BoolScalarParam;
import beast.base.spec.inference.parameter.BoolVectorParam;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.IntVectorParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.pkgmgmt.PackageManager;
import phylonco.beast.evolution.likelihood.ReadCountTreeLikelihood;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Post-run genotype sampling for a {@link ReadCountTreeLikelihood} analysis.
 *
 * <pre>
 *   StochasticMapper -xml analysis.xml -log analysis.log -trees analysis.trees -out PREFIX
 *                    [-burnin 10] [-thin 1] [-draws 1] [-seed 1] [-likelihood ID] [-tolerance 1e-6]
 * </pre>
 *
 * For each retained sample of the trace log (after burnin percent, every thin-th) it sets every state
 * node of the analysis to the logged values and the tree to the logged tree of that sample, recomputes
 * the posterior (and stops if it differs from the logged posterior by more than the tolerance), and
 * draws genotypes for every node and site from their posterior given the reads
 * ({@link ReadCountTreeLikelihood#sampleGenotypes}). It writes
 * <ul>
 * <li>PREFIX.genotypes.log: per draw the sampled genotypes of each cell (one column per cell, as the
 *     SampledGenotypeLogger), and the number of branch-site genotype changes the draw implies;</li>
 * <li>PREFIX.genotype_probs.tsv: per cell and site, the proportion of draws with each genotype.</li>
 * </ul>
 */
public class StochasticMapper {

    public static void main(String[] args) throws Exception {
        Map<String, String> a = new HashMap<>();
        for (int i = 0; i + 1 < args.length; i += 2) {
            if (!args[i].startsWith("-")) usage("unexpected argument " + args[i]);
            a.put(args[i].substring(1), args[i + 1]);
        }
        if (args.length % 2 != 0) usage("missing value after " + args[args.length - 1]);
        for (String k : new String[]{"xml", "log", "trees", "out"}) if (!a.containsKey(k)) usage("missing -" + k);
        new StochasticMapper().run(new File(a.get("xml")), new File(a.get("log")), new File(a.get("trees")),
                a.get("out"), Double.parseDouble(a.getOrDefault("burnin", "10")),
                Integer.parseInt(a.getOrDefault("thin", "1")), Integer.parseInt(a.getOrDefault("draws", "1")),
                Long.parseLong(a.getOrDefault("seed", "1")), a.get("likelihood"),
                Double.parseDouble(a.getOrDefault("tolerance", "1e-6")));
    }

    private static void usage(String problem) {
        System.err.println(problem);
        System.err.println("usage: StochasticMapper -xml X -log L -trees T -out PREFIX [-burnin 10] [-thin 1] "
                + "[-draws 1] [-seed 1] [-likelihood ID] [-tolerance 1e-6]");
        System.exit(1);
    }

    public void run(File xml, File log, File treesFile, String out, double burninPercent, int thin, int draws,
                    long seed, String likelihoodID, double tolerance) throws Exception {
        try {
            PackageManager.loadExternalJars();
        } catch (Throwable e) {
            System.err.println("note: could not load external packages (" + e.getMessage() + ")");
        }
        MCMC mcmc = (MCMC) new XMLParser().parseFile(xml);
        State state = mcmc.startStateInput.get();
        Distribution posterior = mcmc.posteriorInput.get();
        ReadCountTreeLikelihood likelihood = findLikelihood(mcmc, likelihoodID);
        Tree tree = (Tree) likelihood.treeInput.get();
        state.initialise();
        state.setPosterior(posterior);

        // the trace: header and rows by sample number
        List<String> header;
        List<double[]> rows = new ArrayList<>();
        try (BufferedReader r = new BufferedReader(new FileReader(log))) {
            String line;
            header = null;
            while ((line = r.readLine()) != null) {
                if (line.startsWith("#") || line.isBlank()) continue;
                String[] f = line.split("\t");
                if (header == null) {
                    header = List.of(f);
                    continue;
                }
                double[] v = new double[f.length];
                for (int i = 0; i < f.length; i++) v[i] = Double.parseDouble(f[i]);
                rows.add(v);
            }
        }
        if (header == null) throw new IOException(log + ": no header");
        Map<String, Integer> col = new HashMap<>();
        for (int i = 0; i < header.size(); i++) col.put(header.get(i), i);

        // logged trees by sample number
        NexusParser nexus = new NexusParser();
        nexus.parseFile(treesFile);
        Map<Long, Tree> trees = new HashMap<>();
        for (Tree t : nexus.trees) {
            java.util.regex.Matcher m = java.util.regex.Pattern.compile("(\\d+)$").matcher(t.getID());
            if (!m.find()) throw new IOException(treesFile + ": tree " + t.getID() + " has no sample number");
            trees.put(Long.parseLong(m.group(1)), t);
        }

        // log columns of every state node except the tree
        Map<StateNode, String[]> columns = new LinkedHashMap<>();
        for (StateNode sn : state.stateNodeInput.get()) {
            if (sn == tree) continue;
            String[] names = loggedColumns(sn);
            if (names.length == 0 || (names.length == 1 && names[0].isEmpty())) {
                // logs nothing (e.g. a fixed root-sequence alignment): kept at its value in the XML
                System.err.println("StochasticMapper: state node " + sn.getID() + " is not in the log; kept at its "
                        + "value in the XML");
                continue;
            }
            for (String n : names) {
                if (!col.containsKey(n)) {
                    throw new IllegalArgumentException("state node " + sn.getID() + ": column " + n + " not in " + log);
                }
            }
            columns.put(sn, names);
        }
        Integer posteriorCol = col.get("posterior");

        Random random = new Random(seed);
        TreeNodes nodes = new TreeNodes(tree, likelihood);
        int nSites = likelihood.dataInput.get().getPatternCount();
        int nCells = nodes.cells.length;
        int S = likelihood.dataInput.get().getDataType().getStateCount();
        long[][][] counts = new long[nCells][nSites][S];
        int[][] states = new int[nSites][tree.getNodeCount()];
        long nDraws = 0;
        int start = (int) Math.ceil(burninPercent / 100.0 * rows.size());
        try (PrintWriter g = new PrintWriter(out + ".genotypes.log")) {
            g.print("Sample\tdraw\tchanges");
            for (String c : nodes.cells) g.print("\tgenotype(" + c + ")");
            g.println();
            for (int i = start; i < rows.size(); i += thin) {
                double[] row = rows.get(i);
                long sample = (long) row[col.get("Sample")];
                Tree logged = trees.get(sample);
                if (logged == null) throw new IllegalArgumentException("no tree for sample " + sample + " in " + treesFile);
                for (Map.Entry<StateNode, String[]> e : columns.entrySet()) set(e.getKey(), e.getValue(), row, col);
                tree.assignFrom(logged);
                double logP = state.robustlyCalcPosterior(posterior);
                if (posteriorCol != null && Math.abs(logP - row[posteriorCol]) > tolerance * Math.max(1.0, Math.abs(logP))) {
                    throw new IllegalStateException("sample " + sample + ": recomputed posterior " + logP
                            + " differs from the logged " + row[posteriorCol]);
                }
                for (int d = 0; d < draws; d++) {
                    likelihood.sampleGenotypes(random::nextDouble, states);
                    long changes = 0;
                    for (int k = 0; k < nSites; k++) {
                        for (Node n : tree.getNodesAsArray()) {
                            if (!n.isRoot() && states[k][n.getNr()] != states[k][n.getParent().getNr()]) changes++;
                        }
                    }
                    g.print(sample + "\t" + d + "\t" + changes);
                    for (int c = 0; c < nCells; c++) {
                        int nr = nodes.leafNr(tree, c);
                        StringBuilder sb = new StringBuilder(nSites);
                        for (int k = 0; k < nSites; k++) {
                            sb.append(nodes.genotypeChar(states[k][nr]));
                            counts[c][k][states[k][nr]]++;
                        }
                        g.print("\t" + sb);
                    }
                    g.println();
                    nDraws++;
                }
            }
        }
        try (PrintWriter p = new PrintWriter(out + ".genotype_probs.tsv")) {
            p.print("cell\tsite");
            for (int s = 0; s < S; s++) p.print("\t" + nodes.genotypeName(s));
            p.println();
            for (int c = 0; c < nCells; c++) {
                for (int k = 0; k < nSites; k++) {
                    p.print(nodes.cells[c] + "\t" + k);
                    for (int s = 0; s < S; s++) p.printf("\t%.6g", counts[c][k][s] / (double) nDraws);
                    p.println();
                }
            }
        }
        System.err.println("StochasticMapper: " + nDraws + " draws from " + ((rows.size() - start + thin - 1) / thin)
                + " samples; wrote " + out + ".genotypes.log and " + out + ".genotype_probs.tsv");
    }

    /** The ReadCountTreeLikelihood of the analysis (the one with this ID, or the only one). */
    private static ReadCountTreeLikelihood findLikelihood(BEASTInterface root, String id) {
        List<ReadCountTreeLikelihood> found = new ArrayList<>();
        collect(root, found, new IdentityHashMap<>());
        if (id != null) found.removeIf(l -> !id.equals(l.getID()));
        if (found.size() != 1) {
            throw new IllegalArgumentException("expected one ReadCountTreeLikelihood" + (id != null ? " with ID " + id : "")
                    + ", found " + found.size() + (id == null ? " (use -likelihood ID)" : ""));
        }
        return found.get(0);
    }

    private static void collect(BEASTInterface o, List<ReadCountTreeLikelihood> found, Map<BEASTInterface, Boolean> seen) {
        if (o == null || seen.put(o, true) != null) return;
        if (o instanceof ReadCountTreeLikelihood l) found.add(l);
        for (BEASTInterface child : o.listActiveBEASTObjects()) collect(child, found, seen);
    }

    /** The column names the state node writes when logged (its Loggable header). */
    private static String[] loggedColumns(StateNode sn) {
        if (!(sn instanceof Loggable)) throw new IllegalArgumentException(sn.getID() + " is not loggable");
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(bytes);
        ((Loggable) sn).init(ps);
        ps.flush();
        return bytes.toString().trim().split("\\s+");
    }

    private static void set(StateNode sn, String[] names, double[] row, Map<String, Integer> col) {
        if (sn instanceof RealScalarParam<?> p) {
            p.set(row[col.get(names[0])]);
        } else if (sn instanceof RealVectorParam<?> p) {
            for (int i = 0; i < names.length; i++) p.set(i, row[col.get(names[i])]);
        } else if (sn instanceof IntScalarParam<?> p) {
            p.set((int) Math.round(row[col.get(names[0])]));
        } else if (sn instanceof IntVectorParam<?> p) {
            for (int i = 0; i < names.length; i++) p.set(i, (int) Math.round(row[col.get(names[i])]));
        } else if (sn instanceof BoolScalarParam p) {
            p.set(row[col.get(names[0])] != 0);
        } else if (sn instanceof BoolVectorParam p) {
            for (int i = 0; i < names.length; i++) p.set(i, row[col.get(names[i])] != 0);
        } else {
            throw new IllegalArgumentException("state node " + sn.getID() + " of type " + sn.getClass().getSimpleName()
                    + " cannot be set from the log");
        }
    }

    /** Cells (leaf taxa in read-count order) and genotype codes. */
    private static final class TreeNodes {
        final String[] cells;
        final beast.base.evolution.datatype.DataType datatype;

        TreeNodes(Tree tree, ReadCountTreeLikelihood l) {
            String[] names = l.readCountModelInput.get().getReadCount().getTaxonNames();
            cells = new String[names.length];
            for (int i = 0; i < names.length; i++) cells[i] = names[i].trim();
            datatype = l.dataInput.get().getDataType();
        }

        int leafNr(Tree tree, int cell) {
            for (Node leaf : tree.getExternalNodes()) if (leaf.getID().equals(cells[cell])) return leaf.getNr();
            throw new IllegalArgumentException("cell " + cells[cell] + " is not a leaf of the tree");
        }

        char genotypeChar(int state) {
            return datatype.getCharacter(state).charAt(0);
        }

        private static final String[] GT10 = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"};

        String genotypeName(int state) {
            return datatype.getStateCount() == 10 ? GT10[state] : datatype.getCharacter(state);
        }
    }
}
