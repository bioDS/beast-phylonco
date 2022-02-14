// This is an empty umbrella build including all the component builds.
// This build is not necessarily needed. The component builds work independently.

rootProject.name = "beast-phylonco"

include("lphy")
include("beast2")
include("lphybeast")

pluginManagement {
    // the repos to load Gradle plugins
    repositories {
        mavenCentral()
//        maven {
//            // to local build/plugins
//            url = uri("${rootDir.parent}/GradlePlugins/platforms/build/releases/")
//            println("Temp repo : ${url}")
//        }
        // add sonatype snapshots repository
        maven {
            url=uri("https://s01.oss.sonatype.org/content/repositories/snapshots/")
        }
        gradlePluginPortal()
    }
}

// https://docs.gradle.org/current/userguide/build_cache.html
// https://docs.gradle.org/current/userguide/build_cache_use_cases.html
buildCache {
    local {
        directory = File(rootDir, "build-cache")
        removeUnusedEntriesAfterDays = 30
        println("Creating local build cache : ${directory}")
    }
}
