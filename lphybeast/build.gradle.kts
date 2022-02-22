plugins {
    `java-library`
    distribution
    `maven-publish`
    signing
    id("io.github.linguaphylo.platforms.lphy-publish") version "0.1.2"
}

// version has to be manually adjusted to keep same between version.xml and here
//version = "0.0.6-SNAPSHOT"
base.archivesName.set("phylonco-lb")

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17
    withSourcesJar()
}

val zippedConfig by configurations.creating

// if the project dependencies ues impl, then impl(proj(..)) will only have source code,
// which is equivalent to project-version.jar.
// if api is used, then all dependencies will pass to here,
// but do not use api unless you have to.
dependencies {
    //*** phylonco lphy + lphy ***//
    implementation("io.github.linguaphylo:lphy:1.2.0")
//    implementation("io.github.bioDS:?:?")
    implementation(project(":lphy"))

    //*** phylonco beast2 + beast2 + beastlab ***//
    implementation(project(":beast2"))

    //*** lphybeast + ... ***//
    zippedConfig("io.github.linguaphylo:lphybeast:0.3.0")
//    implementation(fileTree("dir" to "${lb.get().outputs.dir("lib")}", "include" to "**/*.jar"))
    implementation(files( { lb.get().extra["lblibs"] } ))

    //TODO add rest of lphybeast dependencies of beast2 part
    // non-modular lphy jar incl. all dependencies
    implementation(fileTree("lib") {
//        include("lphy-*-all.jar")
        // beast2 + beastlab are in subproject beast2
        include("*addon*.jar", "feast-*.jar", "Mascot.*.jar", "SA.*.jar", "SSM.*.jar")
    })

    // tests
    testImplementation("junit:junit:4.13.2")
//    testRuntimeOnly(beast2)
}
tasks.compileJava.get().dependsOn("installLPhyBEAST")

// unzip lphybeast-*.zip to ${buildDir}/lphybeast/
val lb = tasks.register<Sync>("installLPhyBEAST") {
    val outDir = "${buildDir}/lphybeast"
    zippedConfig.resolvedConfiguration.resolvedArtifacts.forEach({
        println(name + " --- " + it.file.name)
        if (it.file.name.endsWith("zip")) {
            from(zipTree(it.file))
            into(outDir)
        }
    })
    extra["lblibs"] = fileTree("dir" to "${outDir}/lib", "include" to "**/*.jar")
}

// launch lphybeast
tasks.register("runLPhyBEAST", JavaExec::class.java) {
    // use classpath
    jvmArgs = listOf("-cp", sourceSets.main.get().runtimeClasspath.asPath)
    println("clspath = ${sourceSets.main.get().runtimeClasspath.asPath}")
    mainClass.set("lphybeast.LPhyBEAST")
    setArgs(listOf("-o", "$rootDir/tmp/gt16ErrMod.xml",
        "${project(":lphy").projectDir}/examples/gt16CoalErrModel.lphy"))
}


tasks.jar {
    manifest {
        // shared attr in the root build
        attributes(
            "Implementation-Title" to "Phylonco LPhyBEAST",
            "Implementation-Vendor" to "?",
        )
    }
}

tasks.getByName<Tar>("distTar").enabled = false
// exclude start scripts
//tasks.getByName<CreateStartScripts>("startScripts").enabled = false

// dist as a beast2 package:
// 1. all released b2 packages are excluded;
// 2. lphy-*-all.jar is excluded, because SPI is not working with BEAST;
// 3. cannot use modular jar, because BEAST uses a customised non-module system.
distributions {
    main {
        distributionBaseName.set(project.base.archivesName.get())
        contents {
//            eachFile {  println(relativePath)  }
            includeEmptyDirs = false
            into("lib") {
                // include beast2 part
                from(project(":beast2").tasks.jar)
                //TODO have to include lphy part, e.g. for lphybeast Unit tests
                from(project(":lphy").tasks.jar)
                // include lphybeast part
                from(tasks.jar)
            }
            into("."){
                from("${projectDir}"){ include("version.xml") }
                from("$rootDir") {
                    include("README.md")
                    include("LICENSE")
                }
            }
            // include src jar
            into("src") {
                from(tasks.getByName<Jar>("sourcesJar"))
                from(project(":beast2").tasks.getByName<Jar>("sourcesJar"))
            }
            into("examples") {
                from("${project(":beast2").projectDir}/examples")
            }
        }
    }
}

// beast 2 will remove version from Zip file name, and then decompress
// rm lphybeast-$version from the relative path of files inside Zip to make it working
tasks.withType<Zip>() {
    doFirst {
        if ( name.equals("distZip") ) {
            // only activate in distZip, otherwise will affect all jars and zips,
            // e.g. main class not found in lphybeast-$version.jar.
            eachFile {
                relativePath = RelativePath(true, *relativePath.segments.drop(1).toTypedArray())
                println(relativePath)
            }
        }
    }
}


publishing {
    publications {
        // project.name contains "lphy" substring
        create<MavenPublication>(project.name) {
            artifactId = project.base.archivesName.get()
            artifact(tasks.distZip.get())
            pom {
                description.set("?")
                packaging = "zip"
                developers {
                    developer {
                        name.set("Kylie Chen")
                    }
                }
            }
        }
    }
}

