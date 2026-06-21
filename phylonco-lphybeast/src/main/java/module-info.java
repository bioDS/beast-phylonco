open module phylonco.lphybeast {
    requires lphy.beast;
    requires phylonco.beast;
    requires phylonco.lphy;
    requires beast.base;
    requires lphy.base;
    requires flc;
    requires beast.pkgmgmt;
    requires com.google.common;
    requires beast.classic;

    exports phylonco.lphybeast.spi;
    exports phylonco.lphybeast.tobeast.generators;

    provides lphybeast.spi.LPhyBEASTMapping with phylonco.lphybeast.spi.LBPhylonco;
}