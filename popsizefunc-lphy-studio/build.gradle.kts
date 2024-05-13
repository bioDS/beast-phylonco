plugins {
    `java-library`
    `maven-publish`
    signing
    id("io.github.linguaphylo.platforms.lphy-java") version "0.1.2"
    id("io.github.linguaphylo.platforms.lphy-publish") version "0.1.2"
}

// overwrite version
version = "0.0.1-SNAPSHOT"//-SNAPSHOT"

dependencies {
    implementation(project(":popsizefunc-lphy"))

    api("io.github.linguaphylo:lphy-studio:1.5.1-SNAPSHOT") //-SNAPSHOT

    testImplementation("org.junit.jupiter:junit-jupiter:5.9.2")
}

val maincls : String = "lphystudio.app.LinguaPhyloStudio"
//application {
//    // equivalent to -m lphystudio
//    // need both mainModule and mainClass
//    mainModule.set("lphystudio")
//    // if only mainClass, it will auto add maincls to the end of CMD
//    mainClass.set(maincls)
//    // applicationDefaultJvmArgs = listOf("-Dgreeting.language=en")
//}

// make studio app locating the correct parent path of examples sub-folder
tasks.withType<JavaExec>() {
    // set version into system property
    systemProperty("lphy.studio.version", version)
    // projectDir = ~/WorkSpace/linguaPhylo/lphy-studio/
    // rootDir = projectDir.parent = ~/WorkSpace/linguaPhylo
    // user.dir = ~/WorkSpace/linguaPhylo/, so examples can be loaded properly
    doFirst {
        // equivalent to: java -p ...
        // user.dir=rootDir (~/WorkSpace/linguaPhylo/), so examples can be loaded properly
        jvmArgs = listOf("-p", classpath.asPath, "-Duser.dir=${rootDir}")
        classpath = files()
    }
    doLast {
        println("JavaExec : $jvmArgs")
    }
}

val developers = "LPhy developer team"
tasks.jar {
    manifest {
        // shared attr in the root build
        attributes(
            "Main-Class" to maincls,
            "Implementation-Title" to "Pop-size function",
            "Implementation-Vendor" to developers,
        )
    }
}

publishing {
    publications {
        // project.name contains "lphy" substring
        create<MavenPublication>(project.name) {
            artifactId = project.base.archivesName.get()
            pom {
                description.set("The LPhy extension including pop-size functions.")
                developers {
                    developer {
                        name.set(developers)
                    }
                }
            }
        }

    }
}


/**
 * For LPhy core, set working directory: ~/WorkSpace/linguaPhylo/lphy/doc,
 * and args[0] = version.
 * For extension, set working directory: ~/WorkSpace/$REPO/lphy/doc,
 * and args[0] = version, args[1] = extension name (no space),
 * args[2] = class name to implement LPhyExtension.
 * e.g. args = 0.0.5 "LPhy Extension Phylonco" phylonco.lphy.spi.Phylonco
 *
 * The docs will output to working dir, "user.dir"
 * This is equivalent to: java -p $LPHY/lib -m lphy/lphy.doc.GenerateDocs 1.1.0
 */
val lphyDoc = tasks.register("lphyDoc", JavaExec::class.java) {
    description = "Create LPhy doc"
    dependsOn("assemble")
//    println("user.dir = " + System.getProperty("user.dir"))

    // equivalent to: java -p ...
    jvmArgs = listOf("-p", sourceSets.main.get().runtimeClasspath.asPath,
        // set output to .../lphy/doc
        "-Duser.dir=${layout.projectDirectory.dir("doc")}")

    // -m lphystudio/lphystudio.app.docgenerator.GenerateDocs
    mainModule.set("lphystudio")
    mainClass.set("lphystudio.app.docgenerator.GenerateDocs")
    // such as 1.1.0
    setArgs(listOf("$version"))
}

// junit tests, https://docs.gradle.org/current/dsl/org.gradle.api.tasks.testing.Test.html
tasks.test {
    useJUnitPlatform() {
        excludeTags("dev")
    }
    // set heap size for the test JVM(s)
    minHeapSize = "128m"
    maxHeapSize = "1G"
    // show standard out and standard error of the test JVM(s) on the console
    testLogging.showStandardStreams = true

    reports {
        junitXml.apply {
            isOutputPerTestCase = true // defaults to false
            mergeReruns.set(true) // defaults to false
        }
    }
}
