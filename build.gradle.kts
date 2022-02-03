
// note kotlin js currently depends on node 14, which doesn't exist for
// arm based mac, so build on 4x4 until kotlin switches to node 16

plugins {
     kotlin("js") version "1.6.10"
}

group = "com.oregondsp"
version = "1.3-SNAPSHOT"

kotlin {
    js {
        browser {
        }
        binaries.executable()
    }
}

repositories {
        mavenCentral()
        maven("https://maven.pkg.jetbrains.space/public/p/kotlinx-html/maven")
    }

/*
build.doLast {
    copy {
        from "${buildDir}/lib/oregondsp.js"
        from "${buildDir}/lib/oregondsp.js.map"
        from "${buildDir}/lib/oregondsp.meta.js"
        into "${projectDir}/lib"
    }
}

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
