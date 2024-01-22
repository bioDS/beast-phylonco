package phylonco.beast;

import beast.pkgmgmt.BEASTClassLoader;
import beast.pkgmgmt.PackageManager;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Set;


/**
 * @author Walter Xie
 */
public class TestUtils {

    private TestUtils() {
    }

    /**
     * user.dir/../version.xml
     */
    public static void loadServices() {
        String wd = System.getProperty("user.dir");
        Path projRoot = Paths.get(wd).toAbsolutePath().getParent();
        File dir = new File(projRoot + "/phylonco-beast/lib");
        if (!dir.exists())
            throw new IllegalArgumentException("Cannot find dir " + dir.getAbsolutePath());

        Path vfPath = Paths.get(projRoot.toString(), "phylonco-lphybeast", "version.xml");
        // do not need jar for local source code
        BEASTClassLoader.addServices(vfPath.toString());
//        File builtJar = new File(projRoot + "/phylonco-beast/build/libs/" + "phylonco-beast-0.0.7.jar");
//        loadPackage(vfPath.toFile(), builtJar);

        // require jar files
        File versionF = new File(dir.getPath() + File.separator + "BEAST.base.version.xml");
        File jarF = new File(dir.getPath() + File.separator + "BEAST.base-2.7.6.jar");
        loadPackage(versionF, jarF);
        versionF = new File(dir.getPath() + File.separator + "BEASTlabs.version.xml");
        jarF = new File(dir.getPath() + File.separator + "BEASTlabs.v2.0.2.jar");
        loadPackage(versionF, jarF);
    }

    private static void loadPackage(File versionF, File jarF) {
        if (!versionF.exists())
            throw new IllegalArgumentException("Cannot find version file: " + versionF.getAbsolutePath());
        if (!jarF.exists())
            throw new IllegalArgumentException("Cannot find jar file: " + jarF.getAbsolutePath());
        // print name and version of package
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        Document doc = null;
        try {
            doc = factory.newDocumentBuilder().parse(versionF);
        } catch (SAXException | IOException | ParserConfigurationException e) {
            throw new RuntimeException(e);
        }
        Element packageElement = doc.getDocumentElement();
        String packageName = packageElement.getAttribute("name");
        Map<String, Set<String>> services = PackageManager.parseServices(doc);

        try {
            PackageManager.addURL(jarF.toURI().toURL(), packageName, services);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static Path getFileForResources(String fileName) {
        System.out.println("WD = " + System.getProperty("user.dir"));
        Path fPath = Paths.get("src", "test", "resources", fileName);
        System.out.println("Input file = " + fPath.toAbsolutePath());
        return fPath;
    }

}
