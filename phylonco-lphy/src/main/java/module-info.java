/**
 * @author Walter Xie
 */
module phylonco.lphy {
    requires transitive lphy.core;
    requires transitive lphy.base;

    exports phylonco.lphy.evolution.alignment;
    exports phylonco.lphy.evolution.datatype;
    exports phylonco.lphy.evolution.substitutionmodel;


    // declare what service interface the provider intends to use
    uses lphy.core.spi.Extension;
    provides lphy.core.spi.Extension with phylonco.lphy.spi.PhyloncoImpl, phylonco.lphy.spi.SequenceTypePhyloncoImpl;
}