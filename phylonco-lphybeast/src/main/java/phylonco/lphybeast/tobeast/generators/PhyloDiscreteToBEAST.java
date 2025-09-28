package phylonco.lphybeast.tobeast.generators;

import NestedBD.evolution.likelihood.DiploidOriginLikelihood;
import beast.base.core.BEASTInterface;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.generators.PhyloCTMCToBEAST;
import phylonco.lphy.evolution.copynumbermodel.MarkovTraitEvolution;
import phylonco.lphy.evolution.copynumbermodel.PhyloDiscrete;
import phylonco.lphybeast.tobeast.values.CopyNumberBDToBEAST;

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
        // Create the DiploidOriginLikelihood
        DiploidOriginLikelihood likelihood = new DiploidOriginLikelihood();

        // Get the beast integer alignment
        if (value instanceof beast.base.evolution.alignment.Alignment data) {
            // TODO: check data is type integer alignment
            likelihood.setInputValue("data", data); // Set alignment data
        } else {
            throw new IllegalArgumentException("Require alignment data");
        }
        // Creates tree and clock rate
        constructTreeAndBranchRate(generator, likelihood, context);

        // Set origin time
        likelihood.setInputValue("origtime", new RealParameter("0.0"));

        // Use the default nstates
        int nstates = CopyNumberBDToBEAST.DEFAULT_NSTATES;
        likelihood.setInputValue("nstates", new RealParameter(String.valueOf(nstates)));

        // Get the model
        Value<MarkovTraitEvolution<Integer>> modelValue = generator.getModel();
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
