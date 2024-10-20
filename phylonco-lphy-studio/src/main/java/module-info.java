import readcountmodel.lphystudio.viewer.Integer2DMatrixViewer;
import readcountmodel.lphystudio.viewer.ReadCountModelViewer;

module popsizefunc.lphy.studio {

    requires transitive lphystudio;
    requires phylonco.lphy;

    // Viewer SPI
    uses lphystudio.app.graphicalmodelpanel.viewer.Viewer;
    // declare what service interface the provider intends to use
    provides lphystudio.app.graphicalmodelpanel.viewer.Viewer with ReadCountModelViewer, Integer2DMatrixViewer;

    // Note: to adapt with the system not using Java module but using class path,
    // they need to be declared inside META-INF/services/lphystudio.app.graphicalmodelpanel.viewer.Viewer as well.

}