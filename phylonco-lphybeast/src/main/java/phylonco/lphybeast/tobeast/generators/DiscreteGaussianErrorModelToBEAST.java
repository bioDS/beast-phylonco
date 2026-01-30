package phylonco.lphybeast.tobeast.generators;

import NestedBD.evolution.likelihood.DiploidOriginLikelihoodWithError;
import NestedBD.evolution.substitutionmodel.BD;
import beast.base.core.BEASTInterface;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.inference.operator.kernel.BactrianRandomWalkOperator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.GraphicalModelNode;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.generators.PhyloCTMCToBEAST;
import phylonco.lphy.evolution.copynumbermodel.*;
import NestedBD.evolution.errormodel.DiscreteGaussianErrorModel;

import java.util.Objects;

public class DiscreteGaussianErrorModelToBEAST
        implements GeneratorToBEAST<phylonco.lphy.evolution.copynumbermodel.DiscreteGaussianErrorModel,
        DiploidOriginLikelihoodWithError> {

    @Override
    public DiploidOriginLikelihoodWithError generatorToBEAST(
            phylonco.lphy.evolution.copynumbermodel.DiscreteGaussianErrorModel generator,
            BEASTInterface value,
            BEASTContext context) {

        // Step 1: Get the error alignment (observed data with measurement error)
        assert value instanceof beast.base.evolution.alignment.Alignment;
        beast.base.evolution.alignment.Alignment errAlignment =
                (beast.base.evolution.alignment.Alignment) value;

        // Step 2: Find the PhyloDiscrete generator (contains tree and CopyNumberBD model)
        PhyloDiscrete phyloDiscrete = findPhyloDiscrete(generator);

        // Step 3: Get nstate from the CopyNumberBD model
        int nstate = getNstateFromModel(phyloDiscrete);

        // Step 4: Remove the original likelihood objects from BEAST context
        removeOriginalLikelihood(generator, phyloDiscrete, context);

        // Step 5: Create the BEAST2 DiscreteGaussianErrorModel
        DiscreteGaussianErrorModel beastErrorModel = createBEASTErrorModel(
                generator, errAlignment, nstate, context);

        // Step 6: Create DiploidOriginLikelihoodWithError
        DiploidOriginLikelihoodWithError likelihoodWithError =
                createLikelihoodWithError(errAlignment, beastErrorModel, phyloDiscrete, nstate, context);

        // Step 7: Add to BEAST logging
        context.addExtraLoggable(likelihoodWithError);

        return likelihoodWithError;
    }

    /**
     * Find the PhyloDiscrete generator by tracing back through the alignment input.
     * The PhyloDiscrete contains the tree and the CopyNumberBD model.
     */
    private PhyloDiscrete findPhyloDiscrete(
            phylonco.lphy.evolution.copynumbermodel.DiscreteGaussianErrorModel generator) {

        for (GraphicalModelNode<?> input : Objects.requireNonNull(generator.getInputs())) {
            if (input instanceof Value && input.value() instanceof IntegerCharacterMatrix) {
                Value<IntegerCharacterMatrix> alignmentValue = (Value<IntegerCharacterMatrix>) input;

                if (alignmentValue.getGenerator() instanceof PhyloDiscrete) {
                    return (PhyloDiscrete) alignmentValue.getGenerator();
                }
            }
        }
        throw new IllegalArgumentException("Cannot find PhyloDiscrete!");
    }

    /**
     * Extract nstate from the CopyNumberBD model.
     * This defines the state space {0, 1, ..., nstate-1}.
     */
    private int getNstateFromModel(PhyloDiscrete phyloDiscrete) {
        Value<MarkovTraitEvolution<Integer>> modelValue = phyloDiscrete.getModel();

        if (modelValue.value() instanceof CopyNumberBD copyNumberBD) {
            return copyNumberBD.getNstate().value();
        } else {
            throw new IllegalArgumentException("Expected CopyNumberBD model but found: " +
                    modelValue.value().getClass().getName());
        }
    }

    /**
     * Create the BEAST2 DiscreteGaussianErrorModel with appropriate parameters.
     */
    private DiscreteGaussianErrorModel createBEASTErrorModel(
            phylonco.lphy.evolution.copynumbermodel.DiscreteGaussianErrorModel generator,
            beast.base.evolution.alignment.Alignment errAlignment,
            int nstate,
            BEASTContext context) {

        DiscreteGaussianErrorModel beastErrorModel = new DiscreteGaussianErrorModel();

        // Set datatype from alignment
        DataType beastDataType = errAlignment.getDataType();
        beastErrorModel.setInputValue("datatype", beastDataType);

        // Convert sigma parameter from LPhy to BEAST
        RealParameter sigmaParam = context.getAsRealParameter(generator.getSigma());
        beastErrorModel.setInputValue("sigma", sigmaParam);
        modifyScaleFactor (sigmaParam, context);

        // Set nstate (number of copy number states)
        beastErrorModel.setInputValue("nstate", new IntegerParameter(String.valueOf(nstate)));

        // Initialize the error model (builds error matrix)
        beastErrorModel.initAndValidate();

        return beastErrorModel;
    }

    private void modifyScaleFactor (RealParameter parameter, BEASTContext context) {
        BactrianRandomWalkOperator operator = new BactrianRandomWalkOperator();
        operator.setInputValue("parameter", parameter);
        operator.setInputValue("weight", context.getOperatorWeight(parameter.getDimension() - 1));
        operator.setInputValue("scaleFactor",0.3);
        operator.initAndValidate();
        operator.setID(parameter.getID() + ".deltaExchange");
        // add operator
        context.addExtraOperator(operator);
        // skip default operator schedule
        context.addSkipOperator(parameter);
    }

    /**
     * Create and configure DiploidOriginLikelihoodWithError.
     * This likelihood calculator integrates over the hidden true copy number states.
     */
    private DiploidOriginLikelihoodWithError createLikelihoodWithError(
            beast.base.evolution.alignment.Alignment errAlignment,
            DiscreteGaussianErrorModel errorModel,
            PhyloDiscrete phyloDiscrete,
            int nstate,
            BEASTContext context) {

        DiploidOriginLikelihoodWithError likelihoodWithError = new DiploidOriginLikelihoodWithError();

        // Set the observed alignment (with measurement error)
        likelihoodWithError.setInputValue("data", errAlignment);

        // Set the error model
        likelihoodWithError.setInputValue("errorModel", errorModel);

        // Get and set the tree and branch rate model
        PhyloCTMCToBEAST.constructTreeAndBranchRate(phyloDiscrete, likelihoodWithError, context, false);

        // Set origin time
        likelihoodWithError.setInputValue("origtime", new RealParameter("0.0"));

        // Set number of states
        likelihoodWithError.setInputValue("nstates", new IntegerParameter(String.valueOf(nstate)));

        // Get the BD substitution model and create site model
        Value<MarkovTraitEvolution<Integer>> modelValue = phyloDiscrete.getModel();
        BD bdModel = (BD) context.getBEASTObject(modelValue);

        SiteModel siteModel = new SiteModel();
        siteModel.setInputValue("substModel", bdModel);
        siteModel.initAndValidate();
        likelihoodWithError.setInputValue("siteModel", siteModel);

        // Initialize the likelihood
        likelihoodWithError.initAndValidate();

        return likelihoodWithError;
    }

    /**
     * Remove the original alignment and likelihood that are replaced by the error model version.
     */
    private void removeOriginalLikelihood(
            phylonco.lphy.evolution.copynumbermodel.DiscreteGaussianErrorModel generator,
            PhyloDiscrete phyloDiscrete,
            BEASTContext context) {

        // Find and remove the original true alignment
        for (GraphicalModelNode<?> input : Objects.requireNonNull(generator.getInputs())) {
            if (input instanceof Value && input.value() instanceof IntegerCharacterMatrix) {
                Value<IntegerCharacterMatrix> origAlignmentInput = (Value<IntegerCharacterMatrix>) input;
                BEASTInterface beastOrigAlignment = context.getBEASTObject(origAlignmentInput);

                if (beastOrigAlignment != null) {
                    context.removeBEASTObject(beastOrigAlignment);
                }
                break;
            }
        }

        // Remove the original likelihood
        BEASTInterface originalLikelihood = context.getBEASTObject(phyloDiscrete);
        if (originalLikelihood != null) {
            context.removeBEASTObject(originalLikelihood);
        }
    }

    @Override
    public Class<phylonco.lphy.evolution.copynumbermodel.DiscreteGaussianErrorModel> getGeneratorClass() {
        return phylonco.lphy.evolution.copynumbermodel.DiscreteGaussianErrorModel.class;
    }

    @Override
    public Class<DiploidOriginLikelihoodWithError> getBEASTClass() {
        return DiploidOriginLikelihoodWithError.class;
    }
}