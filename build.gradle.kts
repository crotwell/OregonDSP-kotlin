
// note kotlin js currently depends on node 14, which doesn't exist for
// arm based mac, so build on 4x4 until kotlin switches to node 16

plugins {
     kotlin("js") version "1.6.10"
}

group = "oregondsp"
version = "1.3.0-beta.1"

kotlin {
    js(IR) {
      moduleName = "oregondsp"
      compilations["main"].packageJson {
        customField("description", "Port of OregonDSP library from java to javascript via kotlin.")
        customField("repository", mapOf("type" to "git",
                                        "url" to "https://github.com/crotwell/OregonDSP-kotlin.git"))
        customField("author", "Philip Crotwell <crotwell@seis.sc.edu>")
        customField("license", "LGPL-3.0")

        customField("keywords", listOf(
          "seismology",
          "fft",
          "timeseries",
          "filter",
          "butterworth",
          "chebyshev",
          "seismogram"
        ))
        customField("bugs", mapOf(
          "url" to "https://github.com/crotwell/OregonDSP-kotlin/issues"
        ))
        customField("homepage", "https://github.com/crotwell/OregonDSP-kotlin")
      }
      browser {
      }
      binaries.executable()
    }
}

repositories {
    mavenCentral()
    maven("https://maven.pkg.jetbrains.space/public/p/kotlinx-html/maven")
}


tasks.register<Sync>("copyJsToLib") {
  from ("README.md")
  from (".npmignore")
  from("${buildDir}/js/packages/oregondsp")
  into("${projectDir}/lib")
  dependsOn("browserDevelopmentWebpack")
}
tasks.get("assemble").dependsOn(tasks.get("copyJsToLib"))
/*
build.doLast {
    configurations.compile.each { File file ->
        copy {
            includeEmptyDirs = false

            from zipTree(file.absolutePath)
            into "${buildDir}/kotlinjs-stdlib"
            include { fileTreeElement ->
                def path = fileTreeElement.path
                path.endsWith(".js") && (path.startsWith("META-INF/resources/") || !path.startsWith("META-INF/"))
            }
        }
    }
}
*/
