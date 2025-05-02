package phylonco.lphybeast.tobeast.generators;

import NestedBD.evolution.likelihood.DiploidOriginLikelihood;
import beast.base.core.BEASTInterface;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.lphy.evolution.copynumbermodel.MarkovTraitEvolution;
import phylonco.lphy.evolution.copynumbermodel.PhyloDiscrete;
import phylonco.lphybeast.tobeast.values.CopyNumberBDToBEAST;


public class PhyloDiscreteToBEAST implements GeneratorToBEAST<PhyloDiscrete, DiploidOriginLikelihood> {

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

        // Get the tree
        Value<TimeTree> treeValue = generator.getTree();
        Tree tree = (Tree) context.getBEASTObject(treeValue);
        likelihood.setInputValue("tree", tree); // Set the tree

        // Set origin time
        likelihood.setInputValue("origtime", new RealParameter("0.5"));

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
