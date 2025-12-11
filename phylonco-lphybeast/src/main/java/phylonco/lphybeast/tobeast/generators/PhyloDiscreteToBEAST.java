package phylonco.lphybeast.tobeast.generators;

import NestedBD.evolution.likelihood.DiploidOriginLikelihood;
import beast.base.core.BEASTInterface;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.generators.PhyloCTMCToBEAST;
import phylonco.lphy.evolution.copynumbermodel.CopyNumberBD;
import phylonco.lphy.evolution.copynumbermodel.MarkovTraitEvolution;
import phylonco.lphy.evolution.copynumbermodel.PhyloDiscrete;

/**
 * BEAST converter for PhyloDiscrete models.
 * <p>
 * Converts PhyloDiscrete generative distributions into BEAST's DiploidOriginLikelihood
 * for Bayesian inference of discrete trait evolution.
 */
public class PhyloDiscreteToBEAST implements GeneratorToBEAST<PhyloDiscrete, DiploidOriginLikelihood> {

    // Creates tree and clock rate by reusing PhyloCTMCToBEAST infrastructure.
    public static void constructTreeAndBranchRate(PhyloDiscrete phyloDiscrete, DiploidOriginLikelihood likelihood, BEASTContext context) {
        PhyloCTMCToBEAST.constructTreeAndBranchRate(phyloDiscrete, likelihood, context, false);
    }

    @Override
    public DiploidOriginLikelihood generatorToBEAST(PhyloDiscrete generator, BEASTInterface value, BEASTContext context) {
        // Check type first
        if (!(value instanceof beast.base.evolution.alignment.Alignment alignment)) {
            throw new IllegalArgumentException("Expected alignment data but got: " + value.getClass().getName());
        }

        DiploidOriginLikelihood likelihood = new DiploidOriginLikelihood();
        likelihood.setInputValue("data", alignment);

        // Creates tree and clock rate
        constructTreeAndBranchRate(generator, likelihood, context);

        // Set origin time
        likelihood.setInputValue("origtime", new RealParameter("0.0"));

        // Get the CopyNumberBD model from the generator
        Value<MarkovTraitEvolution<Integer>> modelValue = generator.getModel();
        CopyNumberBD copyNumberBD = (CopyNumberBD) modelValue.value();

        // Get nstates (will use default if not provided by user)
        int nstates = copyNumberBD.getNstate().value();
        likelihood.setInputValue("nstates", new IntegerParameter(String.valueOf(nstates)));

        // Get the BEAST substitution model
        SubstitutionModel evolutionModel = (SubstitutionModel) context.getBEASTObject(modelValue);

        // Create a site model with the substitution model
        beast.base.evolution.sitemodel.SiteModel siteModel = new beast.base.evolution.sitemodel.SiteModel();
        siteModel.setInputValue("substModel", evolutionModel); // Set substModel
        siteModel.initAndValidate();

        // Set siteModel
        likelihood.setInputValue("siteModel", siteModel);
        likelihood.initAndValidate();

        return likelihood;
    }

    @Override
    public Class<PhyloDiscrete> getGeneratorClass() {
        return PhyloDiscrete.class;
    }

    @Override
    public Class<DiploidOriginLikelihood> getBEASTClass() {
        return DiploidOriginLikelihood.class;
    }
}
