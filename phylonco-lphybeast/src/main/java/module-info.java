open module phylonco.lphybeast {
    requires lphy.beast;
    requires phylonco.beast;
    requires phylonco.lphy;
    requires beast.base;
    requires lphy.base;
    requires flc;

    exports phylonco.lphybeast.spi;
    exports phylonco.lphybeast.tobeast.generators;

    provides lphybeast.spi.LPhyBEASTMapping with phylonco.lphybeast.spi.LBPhylonco;
}