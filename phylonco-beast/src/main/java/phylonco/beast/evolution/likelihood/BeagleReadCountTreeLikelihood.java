package phylonco.beast.evolution.likelihood;

import beagle.*;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beast.base.inference.StateNode;
import beast.base.spec.evolution.branchratemodel.StrictClockModel;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.sitemodel.SiteModel;
import mutablealignment.MutableAlignment;
import phylonco.beast.evolution.datatype.NucleotideDiploid10;
import phylonco.beast.evolution.datatype.NucleotideDiploid16;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BeagleReadCountTreeLikelihood extends TreeLikelihood {

    // This property is a comma-delimited list of resource numbers (0 == CPU) to
    // allocate each BEAGLE instance to. If less than the number of instances then
    // will wrap around.
    // note: to use a different device, say device 2, start beast with
    // java -Dbeagle.resource.order=2 beastfx.app.beast.BeastMCMC
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";
    private static final String SCALING_PROPERTY = "beagle.scaling";
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";
    // Which scheme to use if choice not specified (or 'default' is selected):
    private static final BeagleReadCountTreeLikelihood.PartialsRescalingScheme DEFAULT_RESCALING_SCHEME = BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DYNAMIC;

    private static int instanceCount = 0;
    private static List<Integer> resourceOrder = null;
    private static List<Integer> preferredOrder = null;
    private static List<Integer> requiredOrder = null;
    private static List<String> scalingOrder = null;

    private static final int RESCALE_FREQUENCY = 10000;
    private static final int RESCALE_TIMES = 1;

    boolean m_bUseAmbiguities, m_bUseTipLikelihoods;
    int m_nStateCount;
    int m_nNodeCount;

    private double [] currentCategoryRates;
    //    private double [] storedCurrentCategoryRates;
    private double [] currentFreqs;
    private double [] currentCategoryWeights;

    private int invariantCategory = -1;

    final public Input<LikelihoodReadCountModel> readCountModelInput = new Input<>("readCountModel", "read count model to use for partials");

    final public Input<Boolean> useTipsEmpiricalInput = new Input<>("useTipsEmpirical", "use tip ambiguities from data", false);

    protected LikelihoodReadCountModel readCountModel; // read count model here

    protected boolean useTipsEmpirical;

    // alignment-taxon-index -> ReadCount-row-index
    private int[] alignToRC;

    private StateNode[] rcParams;

    public BeagleReadCountTreeLikelihood() {
        dataInput.setRule(Input.Validate.OPTIONAL);
    }

//    // summed log of the per-(tip,site) max-normalisation factors; added back to logP
//    private double globalLogOffset = 0.0;
//    private double storedGlobalLogOffset = 0.0;

    @Override
    public void initAndValidate() {
        boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));
        if (forceJava) {
            return;
        }
        //globalLogOffset = 0.0;
        initialize();
    }

    private boolean initialize() {
        m_nNodeCount = treeInput.get().getNodeCount();
        m_bUseAmbiguities = m_useAmbiguities.get();
        m_bUseTipLikelihoods = m_useTipLikelihoods.get();
        useTipsEmpirical = useTipsEmpiricalInput.get();

        readCountModel = readCountModelInput.get();

        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }

        m_siteModel = (SiteModel.Base) siteModelInput.get();
        substitutionModel = m_siteModel.substModelInput.get();

// infer genotype datatype from substitution model
        int stateCount = substitutionModel.getStateCount();
        DataType genotypeDataType;
        if (stateCount == 16) {
            genotypeDataType = new NucleotideDiploid16();
        } else if (stateCount == 10) {
            genotypeDataType = new NucleotideDiploid10();
        } else {
            throw new IllegalArgumentException(
                    "BeagleReadCountTreeLikelihood expects GT10 or GT16 substitution model, got stateCount=" + stateCount);
        }

// important: set datatype before buildScaffold / initialize
        readCountModel.setDataType(genotypeDataType);

// build scaffold if XML did not provide data
        if (dataInput.get() == null) {
            dataInput.setValue(buildScaffold(genotypeDataType), this);
        }

        m_siteModel.setDataType(dataInput.get().getDataType());

        alignToRC = readCountModel.getAlignToRCIndex();
        rcParams = new StateNode[]{
                readCountModel.epsilonInput.get(),
                readCountModel.deltaInput.get(),
                readCountModel.tInput.get(),
                readCountModel.vInput.get(),
                readCountModel.sInput.get(),
                readCountModel.w1Input.get(),
                readCountModel.w2Input.get(),
        };

        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) {
            branchRateModel = new StrictClockModel();
        }

        m_branchLengths = new double[m_nNodeCount];
        storedBranchLengths = new double[m_nNodeCount];

        m_nStateCount = dataInput.get().getMaxStateCount();
        patternCount = dataInput.get().getPatternCount();

        eigenCount = 1;//this.branchSubstitutionModel.getEigenCount();

        double[] categoryRates = m_siteModel.getCategoryRates(null);
        // check for invariant rates category
        if (m_siteModel.hasPropInvariantCategory) {
            for (int i = 0; i < categoryRates.length; i++) {
                if (categoryRates[i] == 0) {
                    proportionInvariant = m_siteModel.getRateForCategory(i, null);
                    int patterns = dataInput.get().getPatternCount();
                    calcConstantPatternIndices(patterns, stateCount);
                    invariantCategory = i;

                    double [] tmp = new double [categoryRates.length - 1];
                    for (int k = 0; k < invariantCategory; k++) {
                        tmp[k] = categoryRates[k];
                    }
                    for (int k = invariantCategory + 1; k < categoryRates.length; k++) {
                        tmp[k-1] = categoryRates[k];
                    }
                    categoryRates = tmp;
                    break;
                }
            }
            if (getConstantPattern() != null && getConstantPattern().size() > dataInput.get().getPatternCount()) {
                // if there are many more constant patterns than patterns (each pattern can
                // have a number of constant patters, one for each state) it is less efficient
                // to just calculate the TreeLikelihood for constant sites than optimising
                Log.debug("switch off constant sites optimisiation: calculating through separate TreeLikelihood category (as in the olden days)");
                invariantCategory = -1;
                proportionInvariant = 0;
                setConstantPattern(null);
                categoryRates = m_siteModel.getCategoryRates(null);
            }
        }

        this.categoryCount = m_siteModel.getCategoryCount() - (invariantCategory >= 0 ? 1 : 0);
        tipCount = treeInput.get().getLeafNodeCount();

        internalNodeCount = m_nNodeCount - tipCount;

        int compactPartialsCount = tipCount;
        if (m_bUseAmbiguities) {
            // if we are using ambiguities then we don't use tip partials
            compactPartialsCount = 0;
        }

        // one partials buffer for each tip and two for each internal node (for store restore)
        // changed to two buffers for each node
        partialBufferHelper = new BeagleReadCountTreeLikelihood.BufferIndexHelper(m_nNodeCount, 0);

        // two eigen buffers for each decomposition for store and restore.
        eigenBufferHelper = new BeagleReadCountTreeLikelihood.BufferIndexHelper(eigenCount, 0);

        // two matrices for each node less the root
        matrixBufferHelper = new BeagleReadCountTreeLikelihood.BufferIndexHelper(m_nNodeCount, 0);

        // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
        scaleBufferHelper = new BeagleReadCountTreeLikelihood.BufferIndexHelper(getScaleBufferCount(), 0);

        // Attempt to get the resource order from the System Property
        if (resourceOrder == null) {
            resourceOrder = parseSystemPropertyIntegerArray(RESOURCE_ORDER_PROPERTY);
        }
        if (preferredOrder == null) {
            preferredOrder = parseSystemPropertyIntegerArray(PREFERRED_FLAGS_PROPERTY);
        }
        if (requiredOrder == null) {
            requiredOrder = parseSystemPropertyIntegerArray(REQUIRED_FLAGS_PROPERTY);
        }
        if (scalingOrder == null) {
            scalingOrder = parseSystemPropertyStringArray(SCALING_PROPERTY);
        }

        // first set the rescaling scheme to use from the parser
        rescalingScheme = BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DEFAULT;// = rescalingScheme;
        rescalingScheme = DEFAULT_RESCALING_SCHEME;
        int[] resourceList = null;
        long preferenceFlags = 0;
        long requirementFlags = 0;

        if (scalingOrder.size() > 0) {
            this.rescalingScheme = BeagleReadCountTreeLikelihood.PartialsRescalingScheme.parseFromString(
                    scalingOrder.get(instanceCount % scalingOrder.size()));
        }

        if (resourceOrder.size() > 0) {
            // added the zero on the end so that a CPU is selected if requested resource fails
            resourceList = new int[]{resourceOrder.get(instanceCount % resourceOrder.size()), 0};
            if (resourceList[0] > 0) {
                preferenceFlags |= BeagleFlag.PROCESSOR_GPU.getMask(); // Add preference weight against CPU
            }
        }

        if (preferredOrder.size() > 0) {
            preferenceFlags = preferredOrder.get(instanceCount % preferredOrder.size());
        }

        if (requiredOrder.size() > 0) {
            requirementFlags = requiredOrder.get(instanceCount % requiredOrder.size());
        }

        if (scaling.get().equals(Scaling.always)) {
            this.rescalingScheme = BeagleReadCountTreeLikelihood.PartialsRescalingScheme.ALWAYS;
        }
        if (scaling.get().equals(Scaling.none)) {
            this.rescalingScheme = BeagleReadCountTreeLikelihood.PartialsRescalingScheme.NONE;
        }

        // Define default behaviour here
        if (this.rescalingScheme == BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DEFAULT) {
            //if GPU: the default is^H^Hwas dynamic scaling in BEAST, now NONE
            if (resourceList != null && resourceList[0] > 1) {
                //this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
                this.rescalingScheme = BeagleReadCountTreeLikelihood.PartialsRescalingScheme.NONE;
            } else { // if CPU: just run as fast as possible
                //this.rescalingScheme = PartialsRescalingScheme.NONE;
                // Dynamic should run as fast as none until first underflow
                this.rescalingScheme = BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DYNAMIC;
            }
        }

        if (this.rescalingScheme == BeagleReadCountTreeLikelihood.PartialsRescalingScheme.AUTO) {
            preferenceFlags |= BeagleFlag.SCALING_AUTO.getMask();
            useAutoScaling = true;
        } else {
//                preferenceFlags |= BeagleFlag.SCALING_MANUAL.getMask();
        }
        String r = System.getProperty(RESCALE_FREQUENCY_PROPERTY);
        if (r != null) {
            rescalingFrequency = Integer.parseInt(r);
            if (rescalingFrequency < 1) {
                rescalingFrequency = RESCALE_FREQUENCY;
            }
        }

        if (preferenceFlags == 0 && resourceList == null) { // else determine dataset characteristics
            if (m_nStateCount == 4 && patternCount < 10000) // TODO determine good cut-off
                preferenceFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
        }

        if (substitutionModel.canReturnComplexDiagonalization()) {
            requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
        }

        instanceCount++;

        try {
            beagle = BeagleFactory.loadBeagleInstance(
                    tipCount,
                    partialBufferHelper.getBufferCount(),
                    compactPartialsCount,
                    m_nStateCount,
                    patternCount,
                    eigenBufferHelper.getBufferCount(),            // eigenBufferCount
                    matrixBufferHelper.getBufferCount(),
                    categoryCount,
                    scaleBufferHelper.getBufferCount(), // Always allocate; they may become necessary
                    resourceList,
                    preferenceFlags,
                    requirementFlags
            );
        } catch (Exception e) {
            beagle = null;
        }
        if (beagle == null) {
            return false;
        }

        InstanceDetails instanceDetails = beagle.getDetails();
        ResourceDetails resourceDetails = null;

        if (instanceDetails != null) {
            resourceDetails = BeagleFactory.getResourceDetails(instanceDetails.getResourceNumber());
            if (resourceDetails != null) {
                StringBuilder sb = new StringBuilder("  Using BEAGLE version: " + BeagleInfo.getVersion()
                        + " resource ");
                sb.append(resourceDetails.getNumber()).append(": ");
                sb.append(resourceDetails.getName()).append("\n");
                if (resourceDetails.getDescription() != null) {
                    String[] description = resourceDetails.getDescription().split("\\|");
                    for (String desc : description) {
                        if (desc.trim().length() > 0) {
                            sb.append("    ").append(desc.trim()).append("\n");
                        }
                    }
                }
                sb.append("    with instance flags: ").append(instanceDetails.toString());
                Log.info.println(sb.toString());
            } else {
                Log.warning.println("  Error retrieving BEAGLE resource for instance: " + instanceDetails.toString());
                beagle = null;
                return false;
            }
        } else {
            Log.warning.println("  No external BEAGLE resources available, or resource list/requirements not met, using Java implementation");
            beagle = null;
            return false;
        }
        Log.warning.println("  " + (m_bUseAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
        if (readCountModel != null) {
            Log.warning.println("  Using read count model in tip likelihoods");
        } else if (useTipsEmpirical || m_bUseTipLikelihoods) {
            Log.warning.println(" Using character uncertainty in tip likelihoods");
        }
        Log.warning.println("  With " + patternCount + " unique site patterns.");

        Node[] nodes = treeInput.get().getNodesAsArray();
        for (int i = 0; i < tipCount; i++) {
            int taxon = getTaxonIndex(nodes[i].getID(), dataInput.get());
            if (m_bUseAmbiguities || m_bUseTipLikelihoods || useTipsEmpirical || readCountModel != null) {
                setPartials(beagle, nodes[i], taxon, false);
            } else {
                setStates(beagle, i, taxon);
            }
        }

        if (dataInput.get().isAscertained) {
            ascertainedSitePatterns = true;
        }

        double[] patternWeights = new double[patternCount];
        for (int i = 0; i < patternCount; i++) {
            patternWeights[i] = dataInput.get().getPatternWeight(i);
        }
        beagle.setPatternWeights(patternWeights);

        if (this.rescalingScheme == BeagleReadCountTreeLikelihood.PartialsRescalingScheme.AUTO &&
                resourceDetails != null &&
                (resourceDetails.getFlags() & BeagleFlag.SCALING_AUTO.getMask()) == 0) {
            // If auto scaling in BEAGLE is not supported then do it here
            this.rescalingScheme = BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DYNAMIC;
            Log.warning.println("  Auto rescaling not supported in BEAGLE, using : " + this.rescalingScheme.getText());
        } else {
            Log.warning.println("  Using rescaling scheme : " + this.rescalingScheme.getText());
        }

        if (this.rescalingScheme == BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DYNAMIC) {
            everUnderflowed = false; // If false, BEAST does not rescale until first under-/over-flow.
        }

        updateReadCountModel = true;
        updateSubstitutionModel = true;
        updateSiteModel = true;
        // some subst models (e.g. WAG) never become dirty, so set up subst models right now
        setUpSubstModel();
        // set up error model
        if (readCountModel != null) {
            readCountModel.initialize();
        }
        // set up sitemodel
        beagle.setCategoryRates(categoryRates);
        currentCategoryRates = categoryRates;
        currentFreqs = new double[m_nStateCount];
        currentCategoryWeights = new double[categoryRates.length];

        return true;
    }

    /**
     * Drops the internally-built scaffold from the {@code data} input so a generator does not
     * serialise it into the XML (the integrated model has no genotype alignment — the only data is
     * the read counts). BEAST rebuilds the scaffold from the read counts in {@link #initAndValidate}
     * when it parses and initialises the run.
     */
    public void dropScaffoldInput() {
        dataInput.setValue(null, this);
    }

    /**
     * Builds the constant identity-pattern scaffold the parent {@code TreeLikelihood} requires, from
     * the read-count taxa and site count plus the genotype datatype. The genotype values are dummy
     * (only the taxa / site-count / datatype / pattern structure is used); a {@code MutableAlignment}
     * keeps one pattern per site (no compression), which the per-site read-count mapping relies on.
     */
    private Alignment buildScaffold(DataType genotypeDataType) {
        ReadCount rc = readCountModel.getReadCount();
        String[] taxa = rc.getTaxaNames();
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

        MutableAlignment scaffold = new MutableAlignment(dense);
        scaffold.setID("readCountScaffold");
        scaffold.initAndValidate();
        return scaffold;
    }

    private static List<Integer> parseSystemPropertyIntegerArray(String propertyName) {
        List<Integer> order = new ArrayList<>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    int n = Integer.parseInt(part.trim());
                    order.add(n);
                } catch (NumberFormatException nfe) {
                    Log.warning.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    private static List<String> parseSystemPropertyStringArray(String propertyName) {

        List<String> order = new ArrayList<>();

        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    String s = part.trim();
                    order.add(s);
                } catch (NumberFormatException nfe) {
                    Log.warning.println("Invalid getEigenDecompositionentry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }


    protected int getScaleBufferCount() {
        return internalNodeCount + 1;
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param node    node
     * @param taxon the taxon
     * @param flip whether to flip partial buffer index
     */
    protected final void setPartials(Beagle beagle,
                                     Node node, int taxon, boolean flip) {
        Alignment data = dataInput.get();

        int nrOfStates = data.getDataType().getStateCount();
        int nrOfPatterns = data.getPatternCount();
        double[] partials = new double[nrOfPatterns * nrOfStates * categoryCount];

        int alignIdx = getTaxonIndex(node.getID(), data);
        int rcIdx = alignToRC[alignIdx];

        double[] logRC = new double[nrOfStates];
        int v = 0;
        for (int i = 0; i < patternCount; i++) {
            double[] tipProbabilities;
            if (readCountModel != null) { // read count model partials instead of error model
                int[] reads = readCountModel.getReadCount().getReadCounts(rcIdx, i);
                int coverage = readCountModel.getCoverage(taxon, i);
                double max = Double.NEGATIVE_INFINITY;
                for (int g = 0; g < m_nStateCount; g++) {
                    logRC[g] = readCountModel.logLiklihoodRC(g, reads, coverage, rcIdx);
                    if (logRC[g] > max) max = logRC[g];
                }
                //globalLogOffset += max;
                for (int state = 0; state < m_nStateCount; state++) {
                    //partials[v++] = Math.exp(logRC[state] - max);
                    partials[v++] = Math.exp(logRC[state]);
                }
            } else {
                tipProbabilities = data.getTipLikelihoods(taxon, i);
                if (tipProbabilities != null) {
                    for (int state = 0; state < m_nStateCount; state++) {
                        partials[v++] = tipProbabilities[state];
                    }
                } else {
                    int stateIndex = data.getPattern(taxon, i);
                    boolean[] stateSet = data.getStateSet(stateIndex);
                    for (int state = 0; state < m_nStateCount; state++) {
                        if (stateSet[state]) {
                            partials[v++] = 1.0;
                        } else {
                            partials[v++] = 0.0;
                        }
                    }
                }
            }
        }

        // if there is more than one category then replicate the partials for each
        int n = patternCount * m_nStateCount;
        int k = n;
        for (int i = 1; i < categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }
        if (flip) {
            partialBufferHelper.flipOffset(node.getNr());
        }


        beagle.setPartials(partialBufferHelper.getOffsetIndex(node.getNr()), partials);

    }

    public int getPatternCount() {
        return patternCount;
    }

    void setUpSubstModel() {
        // we are currently assuming a no-category model...
        // TODO More efficient to update only the substitution model that changed, instead of all
        for (int i = 0; i < eigenCount; i++) {
            //EigenDecomposition ed = m_substitutionModel.getEigenDecomposition(i, 0);
            EigenDecomposition ed = substitutionModel.getEigenDecomposition(null);

            eigenBufferHelper.flipOffset(i);

            beagle.setEigenDecomposition(
                    eigenBufferHelper.getOffsetIndex(i),
                    ed.getEigenVectors(),
                    ed.getInverseEigenVectors(),
                    ed.getEigenValues());
        }
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param nodeIndex     nodeIndex
     * @param taxon         the taxon
     */
    protected final void setStates(Beagle beagle,
                                   int nodeIndex, int taxon) {
        Alignment data = dataInput.get();
        int i;

        int[] states = new int[patternCount];

        for (i = 0; i < patternCount; i++) {
            int code = data.getPattern(taxon, i);
            int[] statesForCode = data.getDataType().getStatesForCode(code);
            if (statesForCode.length==1)
                states[i] = statesForCode[0];
            else
                states[i] = code; // Causes ambiguous states to be ignored.
        }

        beagle.setTipStates(nodeIndex, states);
    }

    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
            if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
                throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
    }


//    public void setStates(int tipIndex, int[] states) {
//        System.err.println("BTL:setStates");
//        beagle.setTipStates(tipIndex, states);
//        makeDirty();
//    }
//
//    public void getStates(int tipIndex, int[] states) {
//        System.err.println("BTL:getStates");
//        beagle.getTipStates(tipIndex, states);
//    }


    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    public boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;

        if (readCountModel != null && readCountModel.somethingIsDirty()) { // change to update read count model boolean
            updateReadCountModel = true;
            hasDirt = Tree.IS_DIRTY;
            return true;
        }

        for (StateNode p : rcParams) {
            if (p != null && p.somethingIsDirty()) {
                updateReadCountModel = true;
                hasDirt = Tree.IS_FILTHY;
                return true;
            }
        }




        double[] categoryRates = m_siteModel.getCategoryRates(null);
        if (getConstantPattern() != null) {
            double [] tmp = new double [categoryRates.length - 1];
            for (int k = 0; k < invariantCategory; k++) {
                tmp[k] = categoryRates[k];
            }
            for (int k = invariantCategory + 1; k < categoryRates.length; k++) {
                tmp[k-1] = categoryRates[k];
            }
            categoryRates = tmp;
        }
        for (int i = 0; i < categoryRates.length; i++) {
            if (categoryRates[i] != currentCategoryRates[i]) {
                updateSiteModel = true;
                break;
            }
        }
        //updateSiteModel |= m_siteModel.isDirtyCalculation();

        if (substitutionModel instanceof CalculationNode) {
            updateSubstitutionModel |= ((CalculationNode) substitutionModel).somethingIsDirty();
        }

        if (dataInput.get().somethingIsDirty()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.somethingIsDirty()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.somethingIsDirty()) {
            //m_nHasDirt = Tree.IS_FILTHY;
            return true;
        }

        return treeInput.get().somethingIsDirty();
    }

    /**
     * Stores the additional state other than model components
     */
    @Override
    public void store() {
        //storedGlobalLogOffset = globalLogOffset;
        readCountModel.store(); // read count model store
        partialBufferHelper.storeState();
        eigenBufferHelper.storeState();
        matrixBufferHelper.storeState();

        if (useScaleFactors || useAutoScaling) { // Only store when actually used
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
//            storedRescalingCount = rescalingCount;
        }
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }

    @Override
    public void restore() {
        //globalLogOffset = storedGlobalLogOffset;
        updateSiteModel = true; // this is required to upload the categoryRates to BEAGLE after the restore

        readCountModel.restore(); // read count model restore
        partialBufferHelper.restoreState();
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.restoreState();
            int[] tmp2 = storedScaleBufferIndices;
            storedScaleBufferIndices = scaleBufferIndices;
            scaleBufferIndices = tmp2;
//            rescalingCount = storedRescalingCount;
        }

//        updateRestrictedNodePartials = true;
        super.restore();
        //double[] tmp = m_branchLengths;
        //m_branchLengths = storedBranchLengths;
        //storedBranchLengths = tmp;
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    protected double[] getLeafPartials(Node node) {
        Alignment data = dataInput.get();
        int nrOfStates = data.getDataType().getStateCount();
        int nrOfPatterns = data.getPatternCount();
        double[] partials = new double[nrOfPatterns * nrOfStates];
        int t = getTaxonIndex(node.getID(), data); // taxon index
        int rcIdx = alignToRC[t];
        int i = 0;
        for (int p = 0; p < nrOfPatterns; p++) {
            int state = data.getPattern(t, p);
            double[] tipLikelihoods;
            double[] logRC = new double[nrOfStates];
            if (useTipsEmpirical) {
                tipLikelihoods = data.getTipLikelihoods(t, p);
            } else {
                    tipLikelihoods = new double[nrOfStates];
                    // read count model partials instead of error model
                        int[] reads = readCountModel.getReadCount().getReadCounts(rcIdx, state);
                        int coverage = readCountModel.getCoverage(t, state);
                        double max = Double.NEGATIVE_INFINITY;
                        for (int g = 0; g < m_nStateCount; g++) {
                            logRC[g] = readCountModel.logLiklihoodRC(g, reads, coverage, rcIdx);
                            if (logRC[g] > max) max = logRC[g];
                        }
                        //globalLogOffset += max;
                        for (int v = 0; v < m_nStateCount; v++) {
                            //tipLikelihoods[v] = Math.exp(logRC[v] - max);
                            tipLikelihoods[v] = Math.exp(logRC[v]);
                        }
            }
            for (int s = 0; s < nrOfStates; s++) {
                partials[i] = tipLikelihoods[s];
                i++;
            }
        }
        return partials;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {
        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (matrixUpdateIndices == null) {
            matrixUpdateIndices = new int[eigenCount][m_nNodeCount];
            branchLengths = new double[eigenCount][m_nNodeCount];
            branchUpdateCount = new int[eigenCount];
            scaleBufferIndices = new int[internalNodeCount];
            storedScaleBufferIndices = new int[internalNodeCount];
        }

        if (operations == null) {
            operations = new int[1][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
            operationCount = new int[1];
        }

        recomputeScaleFactors = false;

        if (this.rescalingScheme == BeagleReadCountTreeLikelihood.PartialsRescalingScheme.ALWAYS) {
            useScaleFactors = true;
            recomputeScaleFactors = true;
        } else if (this.rescalingScheme == BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DYNAMIC && everUnderflowed) {
            useScaleFactors = true;
            if (rescalingCountInner < RESCALE_TIMES) {
                recomputeScaleFactors = true;
                hasDirt = Tree.IS_FILTHY;// makeDirty();
//                System.err.println("Recomputing scale factors");
            }

            rescalingCountInner++;
            rescalingCount++;
            if (rescalingCount > RESCALE_FREQUENCY) {
                rescalingCount = 0;
                rescalingCountInner = 0;
            }
        } else if (this.rescalingScheme == BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DELAYED && everUnderflowed) {
            useScaleFactors = true;
            recomputeScaleFactors = true;
            hasDirt = Tree.IS_FILTHY;
            rescalingCount++;
        }

        for (int i = 0; i < eigenCount; i++) {
            branchUpdateCount[i] = 0;
        }
        operationListCount = 0;

        operationCount[0] = 0;

        if (updateReadCountModel) {
            readCountModel.initialize();
        }

        final Node root = treeInput.get().getRoot();
        traverse(root, null, true);

        if (updateReadCountModel) {
            updateReadCountModel = false;
        }

        if (updateSubstitutionModel) {
            setUpSubstModel();
        }

        if (updateSiteModel) {
            double[] categoryRates = m_siteModel.getCategoryRates(null);
            if (getConstantPattern() != null) {
                double [] tmp = new double [categoryRates.length - 1];
                for (int k = 0; k < invariantCategory; k++) {
                    tmp[k] = categoryRates[k];
                }
                for (int k = invariantCategory + 1; k < categoryRates.length; k++) {
                    tmp[k-1] = categoryRates[k];
                }
                categoryRates = tmp;
            }
            for (int i = 0; i < categoryRates.length; i++) {
                if (categoryRates[i] != currentCategoryRates[i]) {
                    beagle.setCategoryRates(categoryRates);
                    i = categoryRates.length;
                }
            }
            currentCategoryRates = categoryRates;
        }

        for (int i = 0; i < eigenCount; i++) {
            if (branchUpdateCount[i] > 0) {
                beagle.updateTransitionMatrices(
                        eigenBufferHelper.getOffsetIndex(i),
                        matrixUpdateIndices[i],
                        null,
                        null,
                        branchLengths[i],
                        branchUpdateCount[i]);
            }
        }

//        if (COUNT_TOTAL_OPERATIONS) {
//            for (int i = 0; i < eigenCount; i++) {
//                totalMatrixUpdateCount += branchUpdateCount[i];
//            }
//
//            for (int i = 0; i <= numRestrictedPartials; i++) {
//                totalOperationCount += operationCount[i];
//            }
//        }

        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {

            beagle.updatePartials(operations[0], operationCount[0], Beagle.NONE);

            int rootIndex = partialBufferHelper.getOffsetIndex(root.getNr());

            double[] categoryWeights = m_siteModel.getCategoryProportions(null);
            if (getConstantPattern() != null) {
                double [] tmp = new double [categoryWeights.length - 1];
                for (int k = 0; k < invariantCategory; k++) {
                    tmp[k] = categoryWeights[k];
                }
                for (int k = invariantCategory + 1; k < categoryWeights.length; k++) {
                    tmp[k-1] = categoryWeights[k];
                }
                categoryWeights = tmp;
            }
            double[] frequencies = substitutionModel.getFrequencies();

            int cumulateScaleBufferIndex = Beagle.NONE;
            if (useScaleFactors) {

                if (recomputeScaleFactors) {
                    scaleBufferHelper.flipOffset(internalNodeCount);
                    cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                    beagle.resetScaleFactors(cumulateScaleBufferIndex);
                    beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, cumulateScaleBufferIndex);
                } else {
                    cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                }
            } else if (useAutoScaling) {
                beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, Beagle.NONE);
            }

            // these could be set only when they change but store/restore would need to be considered

            for (int i = 0; i < categoryWeights.length; i++) {
                if (categoryWeights[i] != currentCategoryWeights[i]) {
                    beagle.setCategoryWeights(0, categoryWeights);
                    i = categoryWeights.length;
                }
            }
            currentCategoryWeights = categoryWeights;
            for (int i = 0; i < frequencies.length; i++) {
                if (frequencies[i] != currentFreqs[i]) {
                    beagle.setStateFrequencies(0, frequencies);
                    i = frequencies.length;
                }
            }
            currentFreqs = frequencies;

            double[] sumLogLikelihoods = new double[1];

            beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                    new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

            logL = sumLogLikelihoods[0];

            if (ascertainedSitePatterns) {
                // Need to correct for ascertainedSitePatterns
                beagle.getSiteLogLikelihoods(patternLogLikelihoods);
                logL = getAscertainmentCorrectedLogLikelihood(dataInput.get(),
                        patternLogLikelihoods, dataInput.get().getWeights(), frequencies);
            } else if (invariantCategory >= 0) {
                beagle.getSiteLogLikelihoods(patternLogLikelihoods);
                int [] patternWeights = dataInput.get().getWeights();
                proportionInvariant = m_siteModel.getProportionInvariant();


                for (int k : getConstantPattern()) {
                    int i = k / m_nStateCount;
                    int j = k % m_nStateCount;
                    patternLogLikelihoods[i] = (Math.log(Math.exp(patternLogLikelihoods[i]) + proportionInvariant * frequencies[j]));
                }

                logL = 0.0;
                for (int i = 0; i < patternCount; i++) {
                    logL += patternLogLikelihoods[i] * patternWeights[i];
                }
            }

            if (Double.isNaN(logL) || Double.isInfinite(logL)) {
                everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt && (rescalingScheme == BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DYNAMIC || rescalingScheme == BeagleReadCountTreeLikelihood.PartialsRescalingScheme.DELAYED)) {
                    // we have had a potential under/over flow so attempt a rescaling
                    useScaleFactors = true;
                    recomputeScaleFactors = true;

                    for (int i = 0; i < eigenCount; i++) {
                        branchUpdateCount[i] = 0;
                    }

                    operationCount[0] = 0;

                    // traverse again but without flipping partials indices as we
                    // just want to overwrite the last attempt. We will flip the
                    // scale buffer indices though as we are recomputing them.
                    traverse(root, null, false);

                    done = false; // Run through do-while loop again
                    firstRescaleAttempt = false; // Only try to rescale once
                } else {
                    // we have already tried a rescale, not rescaling or always rescaling
                    // so just return the likelihood...
                    done = true;
                }
            } else {
                done = true; // No under-/over-flow, then done
            }

        } while (!done);

        // If these are needed...
        //beagle.getSiteLogLikelihoods(patternLogLikelihoods);

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
//        for (int i = 0; i < m_nNodeCount; i++) {
//            updateNode[i] = false;
//        }

        updateSubstitutionModel = false;
        updateSiteModel = false;
        //********************************************************************

        //logP = logL + globalLogOffset;
        logP = logL;
        return logP;
    }

//    protected void getPartials(int number, double[] partials) {
//        int cumulativeBufferIndex = Beagle.NONE;
//        /* No need to rescale partials */
//        beagle.getPartials(partialBufferHelper.getOffsetIndex(number), cumulativeBufferIndex, partials);
//    }

    protected void setPartials(int number, double[] partials) {
        beagle.setPartials(partialBufferHelper.getOffsetIndex(number), partials);
    }

    private double getAscertainmentCorrectedLogLikelihood(Alignment patternList,
                                                          double[] patternLogLikelihoods,
                                                          int[] patternWeights,
                                                          double [] frequencies) {
        if (getConstantPattern() != null) {
            proportionInvariant = m_siteModel.getProportionInvariant();
            for (int k : getConstantPattern()) {
                int i = k / m_nStateCount;
                int j = k % m_nStateCount;
                patternLogLikelihoods[i] = (Math.log(Math.exp(patternLogLikelihoods[i]) + proportionInvariant * frequencies[j]));
            }
        }

        double logL = 0.0;
        double ascertainmentCorrection = patternList.getAscertainmentCorrection(patternLogLikelihoods);
        for (int i = 0; i < patternCount; i++) {
            logL += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
        }
        return logL;
    }

    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @param node           node
     * @param operatorNumber operatorNumber
     * @param flip           flip
     * @return boolean
     */
    private int traverse(Node node, int[] operatorNumber, boolean flip) {
        int nodeNum = node.getNr();

        //Node parent = node.getParent();

        if (operatorNumber != null) {
            operatorNumber[0] = -1;
        }

        // First update the transition probability matrix(ices) for this branch
        int update = (node.isDirty() | hasDirt);
//        if (parent!=null) {
//        	update |= parent.isDirty();
//        }
        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeNum])) {
            m_branchLengths[nodeNum] = branchTime;
            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            if (flip) {
                // first flip the matrixBufferHelper
                matrixBufferHelper.flipOffset(nodeNum);
            }

            // then set which matrix to update
            final int eigenIndex = 0;// = m_substitutionModel.getBranchIndex(node);
            final int updateCount = branchUpdateCount[eigenIndex];
            matrixUpdateIndices[eigenIndex][updateCount] = matrixBufferHelper.getOffsetIndex(nodeNum);

//            if (!m_substitutionModel.canReturnDiagonalization()) {
//            	m_substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate, m_fProbabilities);
//            	int matrixIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
//            	beagle.setTransitionMatrix(matrixIndex, m_fProbabilities, 1);
//            }

            branchLengths[eigenIndex][updateCount] = branchTime;
            branchUpdateCount[eigenIndex]++;

            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            Node child1 = node.getLeft();
            final int[] op1 = {-1};
            final int update1 = traverse(child1, op1, flip);

            Node child2 = node.getRight();
            final int[] op2 = {-1};
            final int update2 = traverse(child2, op2, flip);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                int x = operationCount[operationListCount] * Beagle.OPERATION_TUPLE_SIZE;

                if (flip) {
                    // first flip the partialBufferHelper
                    partialBufferHelper.flipOffset(nodeNum);
                }

                final int[] operations = this.operations[operationListCount];

                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

                if (useScaleFactors) {
                    // get the index of this scaling buffer
                    int n = nodeNum - tipCount;

                    if (recomputeScaleFactors) {
                        // flip the indicator: can take either n or (internalNodeCount + 1) - n
                        scaleBufferHelper.flipOffset(n);

                        // store the index
                        scaleBufferIndices[n] = scaleBufferHelper.getOffsetIndex(n);

                        operations[x + 1] = scaleBufferIndices[n]; // Write new scaleFactor
                        operations[x + 2] = Beagle.NONE;

                    } else {
                        operations[x + 1] = Beagle.NONE;
                        operations[x + 2] = scaleBufferIndices[n]; // Read existing scaleFactor
                    }

                } else {

                    if (useAutoScaling) {
                        scaleBufferIndices[nodeNum - tipCount] = partialBufferHelper.getOffsetIndex(nodeNum);
                    }
                    operations[x + 1] = Beagle.NONE; // Not using scaleFactors
                    operations[x + 2] = Beagle.NONE;
                }

                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNr()); // source node 1
                operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNr()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNr()); // source node 2
                operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNr()); // source matrix 2

                operationCount[operationListCount]++;

                update |= (update1 | update2);

            }
        } else if (node.isLeaf() && updateReadCountModel) {
            if (flip) {
                // first flip the partialBufferHelper
                partialBufferHelper.flipOffset(nodeNum);
            }

            int taxon = getTaxonIndex(node.getID(), dataInput.get());
            setPartials(beagle, node, taxon, false);
        }

        return update;

    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private int eigenCount;
    private int[][] matrixUpdateIndices;
    private double[][] branchLengths;
    private int[] branchUpdateCount;
    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private int[][] operations;
    private int operationListCount;
    private int[] operationCount;

    protected BeagleReadCountTreeLikelihood.BufferIndexHelper partialBufferHelper;
    public BeagleReadCountTreeLikelihood.BufferIndexHelper getPartialBufferHelper() {
        return partialBufferHelper;
    }

    private /*final*/ BeagleReadCountTreeLikelihood.BufferIndexHelper eigenBufferHelper;
    protected BeagleReadCountTreeLikelihood.BufferIndexHelper matrixBufferHelper;
    protected BeagleReadCountTreeLikelihood.BufferIndexHelper scaleBufferHelper;

    protected /*final*/ int tipCount;
    protected /*final*/ int internalNodeCount;
    protected /*final*/ int patternCount;

    private BeagleReadCountTreeLikelihood.PartialsRescalingScheme rescalingScheme = DEFAULT_RESCALING_SCHEME;
    private int rescalingFrequency = RESCALE_FREQUENCY;
    protected boolean useScaleFactors = false;
    private boolean useAutoScaling = false;
    private boolean recomputeScaleFactors = false;
    private boolean everUnderflowed = false;
    private int rescalingCount = 0;
    private int rescalingCountInner = 0;


    /**
     * the number of rate categories
     */
    protected int categoryCount;

    /**
     * an array used to transfer tip partials
     */
    protected double[] tipPartials;

    /**
     * the BEAGLE library instance
     */
    protected Beagle beagle;

    public Beagle getBeagle() {return beagle;}

    /**
     * Flag to specify error model has changed
     */
    protected boolean updateReadCountModel;

    /**
     * Flag to specify that the substitution model has changed
     */
    protected boolean updateSubstitutionModel;
    protected boolean storedUpdateSubstitutionModel;

    /**
     * Flag to specify that the site model has changed
     */
    protected boolean updateSiteModel;
    protected boolean storedUpdateSiteModel;

    /**
     * Flag to specify if site patterns are acertained
     */

    private boolean ascertainedSitePatterns = false;

    public class BufferIndexHelper {
        /**
         * @param maxIndexValue the number of possible input values for the index
         * @param minIndexValue the minimum index value to have the mirrored buffers
         */
        BufferIndexHelper(int maxIndexValue, int minIndexValue) {
            this.maxIndexValue = maxIndexValue;
            this.minIndexValue = minIndexValue;

            offsetCount = maxIndexValue - minIndexValue;
            indexOffsets = new int[offsetCount];
            storedIndexOffsets = new int[offsetCount];
        }

        public int getBufferCount() {
            return 2 * offsetCount + minIndexValue;
        }

        void flipOffset(int i) {
            if (i >= minIndexValue) {
                indexOffsets[i - minIndexValue] = offsetCount - indexOffsets[i - minIndexValue];
            } // else do nothing
        }

        public int getOffsetIndex(int i) {
            if (i < minIndexValue) {
                return i;
            }
            return indexOffsets[i - minIndexValue] + i;
        }

        void getIndices(int[] outIndices) {
            for (int i = 0; i < maxIndexValue; i++) {
                outIndices[i] = getOffsetIndex(i);
            }
        }

        void storeState() {
            System.arraycopy(indexOffsets, 0, storedIndexOffsets, 0, indexOffsets.length);

        }

        void restoreState() {
            int[] tmp = storedIndexOffsets;
            storedIndexOffsets = indexOffsets;
            indexOffsets = tmp;
        }

        private final int maxIndexValue;
        private final int minIndexValue;
        private final int offsetCount;

        private int[] indexOffsets;
        private int[] storedIndexOffsets;

    } // class BufferIndexHelper

    public enum PartialsRescalingScheme {
        DEFAULT("default"), // whatever our current favourite default is
        NONE("none"),       // no scaling
        DYNAMIC("dynamic"), // rescale when needed and reuse scaling factors
        ALWAYS("always"),   // rescale every node, every site, every time - slow but safe
        DELAYED("delayed"), // postpone until first underflow then switch to 'always'
        AUTO("auto");       // BEAGLE automatic scaling - currently playing it safe with 'always'
//        KICK_ASS("kickAss"),// should be good, probably still to be discovered

        PartialsRescalingScheme(String text) {
            this.text = text;
        }

        public String getText() {
            return text;
        }

        private final String text;

        public static BeagleReadCountTreeLikelihood.PartialsRescalingScheme parseFromString(String text) {
            for (BeagleReadCountTreeLikelihood.PartialsRescalingScheme scheme : BeagleReadCountTreeLikelihood.PartialsRescalingScheme.values()) {
                if (scheme.getText().compareToIgnoreCase(text) == 0)
                    return scheme;
            }
            return DEFAULT;
        }
    }

    @Override
    public double [] getPatternLogLikelihoods() {
        beagle.getSiteLogLikelihoods(patternLogLikelihoods);
        return patternLogLikelihoods.clone();
    }

}
