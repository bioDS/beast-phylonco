plugins {
    `java-library`
}

// version has to be manually adjusted to keep same between version.xml and here
version = "1.0.2"//-SNAPSHOT

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17
    withSourcesJar()
}

dependencies {
    // api can pass beast jars to lphybeast
    api(fileTree("lib") {
        // beast 2
//        include("*.jar")
        exclude("**/*-src.jar")
        exclude("**/BEAST-app-*.jar")
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
            "Implementation-Vendor" to "Phylonco development team",
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
    // BeagleTreeLikelihoodWithErrorTest uses Beagle
    exclude("**/BeagleTreeLikelihoodWithErrorTest*")
}

