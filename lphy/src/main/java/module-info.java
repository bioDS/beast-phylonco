/**
 * @author Walter Xie
 */
module phylonco.lphy {
    requires transitive lphy;

    exports phylonco.lphy.evolution.alignment;
    exports phylonco.lphy.evolution.datatype;
    exports phylonco.lphy.evolution.substitutionmodel;


    // declare what service interface the provider intends to use
    provides lphy.spi.LPhyExtension with phylonco.lphy.spi.Phylonco;
}