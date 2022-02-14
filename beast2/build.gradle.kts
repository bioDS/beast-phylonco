plugins {
    `java-library`
}

// version has to be manually adjusted to keep same between version.xml and here
//version = "0.0.6-SNAPSHOT"
base.archivesName.set("phylonco-b2")

java {
    sourceCompatibility = JavaVersion.VERSION_16
    targetCompatibility = JavaVersion.VERSION_16
    withSourcesJar()
}

dependencies {
    // api can pass beast jars to lphybeast
    api(fileTree("lib") {
        // beast 2
        include("beast-*.jar")
        include("BEASTlabs-*.jar")
    })

    // tests
    testImplementation("junit:junit:4.13.2")

//    testRuntimeOnly(beast2)
}


tasks.jar {
    manifest {
        // shared attr in the root build
        attributes(
            "Implementation-Title" to "Phylonco BEAST2",
            "Implementation-Vendor" to "?",
        )
    }
}

tasks.test {
    useJUnit()
    // useJUnitPlatform()
    // set heap size for the test JVM(s)
    minHeapSize = "128m"
    maxHeapSize = "1G"
    // show standard out and standard error of the test JVM(s) on the console
    testLogging.showStandardStreams = true
    //testLogging.exceptionFormat = org.gradle.api.tasks.testing.logging.TestExceptionFormat.FULL
}

