plugins {
    `java-library`
    `maven-publish`
    signing
    id("io.github.linguaphylo.platforms.lphy-java") version "0.1.1"
    id("io.github.linguaphylo.platforms.lphy-publish") version "0.1.1"
}

version = "0.0.1-SNAPSHOT"
base.archivesName.set("phylonco-lphy")

dependencies {
    // io.github.linguaphylo
    implementation("io.github.linguaphylo:lphy:1.1.0")

    // launch studio from its jar, but not depend on it
    runtimeOnly("io.github.linguaphylo:lphy-studio:1.1.0")

    testImplementation("junit:junit:4.13.2")
//    testRuntimeOnly("org.junit.jupiter:junit-jupiter-engine:4.13")
}

// launch lphy studio from io.github.linguaphylo:lphy-studio:version
tasks.register("runLPhyStudio", JavaExec::class.java) {
    // use module
    jvmArgs = listOf("-p", sourceSets.main.get().runtimeClasspath.asPath)
    mainModule.set("lphystudio")
    mainClass.set("lphystudio.app.LinguaPhyloStudio")
}

// lphy-$version.jar
tasks.jar {
    manifest {
        // shared attr in the root build
        attributes(
            "Implementation-Title" to "Phylonco",
            "Implementation-Vendor" to "?",
        )
    }
}

publishing {
    publications {
        // project.name contains "lphy" substring
        create<MavenPublication>(project.name) {
            artifactId = project.base.archivesName.get()
            from(components["java"])
            pom {
                packaging = "jar"
                description.set("?")
                developers {
                    developer {
                        name.set("Kylie Chen")
                    }
                }
            }
        }
    }
}

// junit tests, https://docs.gradle.org/current/dsl/org.gradle.api.tasks.testing.Test.html
//tasks.test {
//    useJUnit()
//    // useJUnitPlatform()
//    // set heap size for the test JVM(s)
//    minHeapSize = "128m"
//    maxHeapSize = "1G"
//    // show standard out and standard error of the test JVM(s) on the console
//    testLogging.showStandardStreams = true
//    //testLogging.exceptionFormat = org.gradle.api.tasks.testing.logging.TestExceptionFormat.FULL
//}

/**
 * For LPhy core, set working directory: ~/WorkSpace/linguaPhylo/lphy/doc,
 * and args[0] = version.
 * For extension, set working directory: ~/WorkSpace/beast-phylonco/lphy/doc,
 * and args[0] = version, args[1] = extension name (no space),
 * args[2] = class name to implement LPhyExtension.
 * e.g. args = 0.0.5 "LPhy Extension Phylonco" phylonco.lphy.spi.Phylonco
 *
 * The docs will output to working dir, "user.dir"
 * This is equivalent to: java -p $LPHY/lib -m lphy/lphy.doc.GenerateDocs 1.1.0
 */
//TODO working in lphy-1.1.1
val lphyDoc by tasks.register("lphyDoc", JavaExec::class.java) {
    description = "Create LPhy doc"
    dependsOn("assemble")
    println("user.dir = " + System.getProperty("user.dir"))

    // equivalent to: java -p ...
    jvmArgs = listOf("-p", sourceSets.main.get().runtimeClasspath.asPath,
        // set output to .../lphy/doc
        "-Duser.dir=${layout.projectDirectory.dir("doc")}")

    // -m lphy/lphy.doc.GenerateDocs
    mainModule.set("lphy")
    mainClass.set("lphy.doc.GenerateDocs")
    // such as 1.1.0
    setArgs(listOf("$version", "Phylonco", "phylonco.lphy.spi.Phylonco"))
}



