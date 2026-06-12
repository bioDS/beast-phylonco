module phylonco.beast {
    requires beast.base;
    requires beast.pkgmgmt;

    requires beagle;
    requires org.apache.commons.statistics.distribution;
    requires org.apache.commons.numbers.gamma;
    requires org.apache.commons.math4.core;
    requires org.apache.commons.math4.legacy;

    requires static beast.fx;
    requires static javafx.controls;
    requires java.xml;

    requires jdk.jfr;
    requires MutableAlignment;

    exports phylonco.beast.evolution.datatype;
    exports phylonco.beast.evolution.errormodel;
    exports phylonco.beast.evolution.likelihood;
    exports phylonco.beast.evolution.readcountmodel;
    exports phylonco.beast.evolution.substitutionmodel;

    provides beast.base.core.BEASTInterface with
            phylonco.beast.evolution.datatype.NucleotideDiploid10,
            phylonco.beast.evolution.datatype.NucleotideDiploid16,
            phylonco.beast.evolution.datatype.NucleotideMethylation,
            phylonco.beast.evolution.datatype.Ternary,
            phylonco.beast.evolution.errormodel.BinaryErrorModel,
            phylonco.beast.evolution.errormodel.ErrorModelBase,
            phylonco.beast.evolution.errormodel.GT16ErrorModel,
            phylonco.beast.evolution.likelihood.BeagleTreeLikelihoodWithError,
            phylonco.beast.evolution.likelihood.TreeLikelihoodWithError,
            phylonco.beast.evolution.likelihood.TreeLikelihoodWithErrorFast,
            phylonco.beast.evolution.likelihood.TreeLikelihoodWithErrorSlow,
            phylonco.beast.evolution.likelihood.ReadCountTreeLikelihood,
            phylonco.beast.evolution.likelihood.SampledGenotypeLogger,
            phylonco.beast.evolution.substitutionmodel.BinarySubstitutionModel,
            phylonco.beast.evolution.substitutionmodel.GT16,
            phylonco.beast.evolution.substitutionmodel.GT10,
            phylonco.beast.evolution.substitutionmodel.MethylationHKY,
            phylonco.beast.evolution.substitutionmodel.SiFit2,
            phylonco.beast.evolution.substitutionmodel.SiFit3,
            phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel,
            phylonco.beast.evolution.readcountmodel.GibbsSequenceOperator,
            phylonco.beast.evolution.readcountmodel.GibbsAlignmentOperator,
            phylonco.beast.evolution.datatype.ReadCount;
}