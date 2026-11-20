package phylonco.beast.evolution.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import mutablealignment.MutableAlignment;
import phylonco.beast.evolution.alignment.ReadCountAlignment;
import phylonco.beast.evolution.datatype.NucleotideDiploid10;
import phylonco.beast.evolution.datatype.NucleotideDiploid16;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.beast.evolution.datatype.ReadCountMatrix;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;

import java.util.List;

/**
 * Tree likelihood that integrates the latent single-cell genotype alignment OUT of the
 * posterior analytically, instead of sampling it as a {@code MutableAlignment} StateNode.
 *
 * <p>Each tip's Felsenstein partial-likelihood vector is set to
 * {@code P(read counts at cell t, site j | genotype g)} for every genotype {@code g}, taken
 * from {@link LikelihoodReadCountModel#logLiklihoodRC}. Standard pruning over the GT16/GT10
 * substitution model then marginalises the genotypes. This collapses the former two-factor
 * posterior (MATreeLikelihood over sampled genotypes + LikelihoodReadCountModel) into a single
 * factor and removes the genotype StateNode and all its Gibbs operators.</p>
 *
 * <p>The read-count likelihoods are tiny, so each tip's partial vector is rescaled per-site by
 * its maximum before upload (max-normalisation); the summed log of those per-(tip,site) maxima
 * is held in {@link #globalLogOffset} and added back to {@code logP}.</p>
 *
 * <p>There is no genotype alignment in the model — the only data is the {@code ReadCount}. The
 * {@code data} Alignment that the parent {@link TreeLikelihoodWithErrorFast} needs for taxa / site
 * count / pattern indexing is built internally as an identity-pattern scaffold from the read counts
 * (taxa, site count) and the genotype datatype implied by the GT10/GT16 substitution model; its
 * genotype values are ignored. The read-count model referenced via {@link #readCountModelInput}
 * must NOT also be added to the posterior CompoundDistribution (it would be double-counted); it
 * stays in the model graph as this likelihood's input, so BEAST still drives its store/restore.</p>
 */
@Description("Tree likelihood that integrates single-cell genotypes out via read-count tip partials")
public class ReadCountTreeLikelihood extends TreeLikelihoodWithErrorFast {

    final public Input<LikelihoodReadCountModel> readCountModelInput = new Input<>(
            "readCountModel",
            "read-count model supplying per-genotype likelihoods P(reads | genotype) for the tip partials",
            Input.Validate.REQUIRED);

    private LikelihoodReadCountModel readCountModel;

    public ReadCountTreeLikelihood() {
        // no genotype alignment in the integrated model; the scaffold is built from the read counts
        dataInput.setRule(Input.Validate.OPTIONAL);
    }

    // alignment-taxon-index -> ReadCount-row-index
    private int[] alignToRC;

    // the 7 read-count parameters; a change in any forces a tip-partial rebuild.
    // Typed as StateNode because the beast3 spec params are a mix of RealScalarParam and
    // RealVectorParam (both extend StateNode, which carries somethingIsDirty()).
    private StateNode[] rcParams;

    // summed log of the per-(tip,site) max-normalisation factors; added back to logP
    private double globalLogOffset = 0.0;
    private double storedGlobalLogOffset = 0.0;


    @Override
    public void initAndValidate() {
        readCountModel = readCountModelInput.get();

        // genotype datatype implied by the GT10/GT16 substitution model (10 or 16 states)
        SubstitutionModel substModel =
                ((SiteModelInterface.Base) siteModelInput.get()).getSubstitutionModel();
        int stateCount = substModel.getStateCount();
        if (stateCount != 10 && stateCount != 16) {
            throw new IllegalArgumentException(
                    "ReadCountTreeLikelihood expects a 10- or 16-state genotype substitution model "
                    + "(GT10 / GT16), got " + stateCount);
        }
        DataType genotypeDataType = stateCount == 16 ? new NucleotideDiploid16() : new NucleotideDiploid10();

        // tell the read-count helper the genotype datatype (it has no genotype alignment) and refresh
        // its cached Gamma tables before the first getLeafPartials traversal super.initAndValidate triggers
        readCountModel.setDataType(genotypeDataType);
        readCountModel.initialize();
        alignToRC = readCountModel.getAlignToRCIndex();
        rcParams = new StateNode[]{
                (StateNode) readCountModel.epsilonInput.get(),
                (StateNode) readCountModel.deltaInput.get(),
                (StateNode) readCountModel.tInput.get(),
                (StateNode) readCountModel.vInput.get(),
                (StateNode) readCountModel.sInput.get(),
                (StateNode) readCountModel.w1Input.get(),
                (StateNode) readCountModel.w2Input.get(),
        };

        // build the identity-pattern scaffold the parent TreeLikelihood needs, from the read counts
        // (there is no genotype alignment in the model). Kept internal — not an XML input.
        if (dataInput.get() == null) {
            dataInput.setValue(buildScaffold(genotypeDataType), this);
        } else if (dataInput.get() instanceof ReadCountAlignment) {
            // the substitution model is authoritative for the genotype state count. A
            // ReadCountAlignment declares a datatype of its own (it must, to be an Alignment), so
            // reconcile rather than reject: BEAUti lets the user switch GT16 to GT10 in the Site
            // Model tab after the partition already exists.
            ((ReadCountAlignment) dataInput.get()).setGenotypeDataType(genotypeDataType);
        }

        // super.initAndValidate() builds the initial tip partials via setPartials -> getLeafPartials,
        // accumulating globalLogOffset from 0.
        globalLogOffset = 0.0;
        super.initAndValidate();

        Alignment data = dataInput.get();
        if (data.getPatternCount() != data.getSiteCount()) {
            throw new IllegalArgumentException(
                    "ReadCountTreeLikelihood requires an identity-pattern scaffold (patternCount == "
                    + "siteCount); pattern compression would break the per-site read-count mapping.");
        }
        if (m_siteModel.getProportionInvariant() != 0.0) {
            throw new IllegalArgumentException(
                    "ReadCountTreeLikelihood does not support proportionInvariant > 0 "
                    + "(it would add an unscaled term to constant-pattern root slots). Set it to 0.");
        }
    }


    /**
     * Builds the constant identity-pattern scaffold the parent {@code TreeLikelihood} requires, from
     * the read-count taxa and site count plus the genotype datatype. The genotype values are dummy
     * (only the taxa / site-count / datatype / pattern structure is used); a {@code MutableAlignment}
     * keeps one pattern per site (no compression), which the per-site read-count mapping relies on.
     */
    private Alignment buildScaffold(DataType genotypeDataType) {
        ReadCountMatrix rc = readCountModel.getReadCount();
        String[] taxa = rc.getTaxonNames();
        int nSites = rc.getSiteNumber();
        String genoChar = genotypeDataType.getCharacter(0); // any valid genotype char; value ignored
        StringBuilder sb = new StringBuilder(nSites);
        for (int j = 0; j < nSites; j++) sb.append(genoChar);
        String seqStr = sb.toString();

        Alignment dense = new Alignment();
        for (String taxon : taxa) {
            dense.setInputValue("sequence", new Sequence(taxon.trim(), seqStr));
        }
        dense.setInputValue("userDataType", genotypeDataType);
        dense.initAndValidate();
        TaxonSet taxonSet = new TaxonSet();
        taxonSet.initByName("alignment", dense);

        MutableAlignment scaffold = new MutableAlignment(dense);
        scaffold.setInputValue("taxa", taxonSet);
        scaffold.setID("readCountScaffold");
        scaffold.initAndValidate();
        return scaffold;
    }

    /**
     * Tip partial vector for one leaf: for each site (== pattern, identity), the per-genotype
     * read-count likelihood, max-normalised per site. Accumulates the per-site log-max into
     * {@link #globalLogOffset}.
     */
    @Override
    protected double[] getLeafPartials(Node node) {
        Alignment data = dataInput.get();
        int nrOfStates = data.getDataType().getStateCount();
        int nrOfPatterns = data.getPatternCount();
        double[] partials = new double[nrOfPatterns * nrOfStates];

        int alignIdx = getTaxonIndex(node.getID(), data);
        int rcIdx = alignToRC[alignIdx];

        double[] logRC = new double[nrOfStates];
        int i = 0;
        for (int p = 0; p < nrOfPatterns; p++) {
            int[] reads = readCountModel.getReadCount().getReadCounts(rcIdx, p);
            int coverage = readCountModel.getCoverage(alignIdx, p);
            double max = Double.NEGATIVE_INFINITY;
            for (int g = 0; g < nrOfStates; g++) {
                logRC[g] = readCountModel.logLiklihoodRC(g, reads, coverage, rcIdx);
                if (logRC[g] > max) max = logRC[g];
            }
            globalLogOffset += max;
            for (int g = 0; g < nrOfStates; g++) {
                partials[i++] = Math.exp(logRC[g] - max);
            }
        }
        return partials;
    }

    /**
     * Rebuild and re-upload every leaf's partials from the current read-count parameters,
     * recomputing globalLogOffset from scratch.
     */
    @Override
    public void updateLeafPartials() {
        globalLogOffset = 0.0;
        readCountModel.initialize();
        List<Node> leaves = treeInput.get().getExternalNodes();
        for (Node node : leaves) {
            int nodeId = node.getNr();
            likelihoodCore.setNodePartialsForUpdate(nodeId);
            BeerLikelihoodCore beer = (BeerLikelihoodCore) likelihoodCore;
            beer.setCurrentNodePartials(nodeId, getLeafPartials(node));
        }
    }

    @Override
    public double calculateLogP() {
        // super (Fast) re-runs updateLeafPartials() when the leaf-partial dirty flag is set,
        // which resets and re-accumulates globalLogOffset; otherwise globalLogOffset is unchanged.
        double core = super.calculateLogP();
        logP = core + globalLogOffset;
        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        // any read-count parameter change invalidates the tip partials -> full rebuild
        for (StateNode p : rcParams) {
            if (p != null && p.somethingIsDirty()) {
                updateLeafPartials = true;
                hasDirt = Tree.IS_FILTHY;
                return true;
            }
        }
        // otherwise defer to the inherited tree / siteModel / branchRate / data checks
        return super.requiresRecalculation();
    }

    @Override
    public void store() {
        storedGlobalLogOffset = globalLogOffset;
        super.store();
    }

    @Override
    public void restore() {
        globalLogOffset = storedGlobalLogOffset;
        super.restore();
    }
}
