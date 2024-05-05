import java.util.Properties

import Settings._
import sbtbuildinfo.BuildInfoPlugin.autoImport._

val appProperties = settingKey[Properties]("The application properties")

inThisBuild(
  Seq(
    organization     := "de.proteinevolution",
    organizationName := "Dept. of Protein Evolution, Max Planck Institute for Biology",
    startYear        := Some(2018),
    licenses += ("Apache-2.0", url("https://www.apache.org/licenses/LICENSE-2.0")),
    scalaVersion := "2.13.14"
  )
)

lazy val buildInfoSettings = Seq(
  buildInfoKeys := Seq[BuildInfoKey](
    name,
    version,
    scalaVersion,
    sbtVersion,
    "playVersion" -> play.core.PlayVersion.current
  ),
  buildInfoPackage := "build"
)

lazy val coreSettings = Settings.compileSettings ++ Release.settings

lazy val common = (project in file("modules/common"))
  .settings(
    name := "common",
    libraryDependencies ++= Dependencies.commonDeps,
    Settings.compileSettings,
    disableDocs
  )
  .settings(scalacOptions += "-Ymacro-annotations")
  .dependsOn(base, tel)

lazy val jobs = (project in file("modules/jobs"))
  .commonSettings("de.proteinevolution.jobs")
  .enablePlugins(PlayScala)
  .dependsOn(common, auth, base, clusterApi, tel, tools, ui, util, user, statistics)
  .settings(scalacOptions += "-Ymacro-annotations")
  .settings(addCompilerPlugin(("org.typelevel" %% "kind-projector" % "0.13.0").cross(CrossVersion.full)))
  .disablePlugins(PlayLayoutPlugin)

lazy val user = (project in file("modules/user"))
  .commonSettings("de.proteinevolution.user")
  .enablePlugins(PlayScala)
  .dependsOn(common, base)
  .disablePlugins(PlayLayoutPlugin)

lazy val auth = (project in file("modules/auth"))
  .commonSettings("de.proteinevolution.auth")
  .enablePlugins(PlayScala)
  .dependsOn(common, base, tel, util, user)
  .disablePlugins(PlayLayoutPlugin)

lazy val base = (project in file("modules/base"))
  .commonSettings("de.proteinevolution.base")
  .enablePlugins(PlayScala)
  .disablePlugins(PlayLayoutPlugin)

lazy val cluster = (project in file("modules/cluster"))
  .commonSettings("de.proteinevolution.cluster")
  .enablePlugins(PlayScala)
  .dependsOn(common, base, jobs, clusterApi)
  .disablePlugins(PlayLayoutPlugin)

lazy val clusterApi = (project in file("modules/cluster-api"))
  .commonSettings("de.proteinevolution.cluster.api")
  .enablePlugins(PlayScala)
  .dependsOn(common, base)
  .disablePlugins(PlayLayoutPlugin)

lazy val backend = (project in file("modules/backend"))
  .commonSettings("de.proteinevolution.backend")
  .enablePlugins(PlayScala)
  .dependsOn(common, base, auth, jobs, message, tel, tools, user, statistics)
  .disablePlugins(PlayLayoutPlugin)

lazy val statistics = (project in file("modules/statistics"))
  .commonSettings("de.proteinevolution.statistics")
  .enablePlugins(PlayScala)
  .dependsOn(common, base)
  .disablePlugins(PlayLayoutPlugin)

lazy val ui = (project in file("modules/ui"))
  .commonSettings("de.proteinevolution.ui")
  .enablePlugins(PlayScala, BuildInfoPlugin)
  .dependsOn(common, base, tools)
  .settings(buildInfoSettings)
  .disablePlugins(PlayLayoutPlugin)

lazy val message = (project in file("modules/message"))
  .commonSettings("de.proteinevolution.message")
  .enablePlugins(PlayScala)
  .dependsOn(common, base, auth, cluster, jobs, tools)
  .disablePlugins(PlayLayoutPlugin)

lazy val migrations = (project in file("modules/migrations"))
  .commonSettings("de.proteinevolution.migrations")
  .enablePlugins(PlayScala)
  .settings(scalacOptions -= "-Wunused:imports")
  .disablePlugins(PlayLayoutPlugin)

lazy val tel = (project in file("modules/tel"))
  .commonSettings("de.proteinevolution.tel")
  .enablePlugins(PlayScala)
  .disablePlugins(PlayLayoutPlugin)

lazy val tools = (project in file("modules/tools"))
  .commonSettings("de.proteinevolution.tools")
  .enablePlugins(PlayScala)
  .dependsOn(common)
  .settings(scalacOptions += "-Ymacro-annotations")
  .disablePlugins(PlayLayoutPlugin)

lazy val util = (project in file("modules/util"))
  .commonSettings("de.proteinevolution.util")
  .enablePlugins(PlayScala)
  .dependsOn(common, user)
  .disablePlugins(PlayLayoutPlugin)

lazy val root = (project in file("."))
  .enablePlugins(PlayScala, PlayAkkaHttp2Support, JavaAppPackaging, AutomateHeaderPlugin)
  .dependsOn(
    jobs,
    cluster,
    backend,
    ui,
    message,
    migrations,
    user
  )
  .settings(
    coreSettings,
    name := "mpi-toolkit",
    libraryDependencies ++= (Dependencies.commonDeps ++ Dependencies.testDeps),
    pipelineStages := Seq(digest, gzip)
  )

resolvers += "scalaz-bintray".at("https://dl.bintray.com/scalaz/releases")
resolvers ++= Resolver.sonatypeOssRepos("releases")
resolvers ++= Resolver.sonatypeOssRepos("snapshots")

fork := true // required for "sbt run" to pick up javaOptions
javaOptions += "-Dplay.editor=http://localhost:63342/api/file/?file=%s&line=%s"
Test / fork     := true
Test / logLevel := Level.Info

Test / scalacOptions ++= Seq("-Yrangepos")
Test / testOptions += Tests.Argument(TestFrameworks.ScalaTest, "-oD")

PlayKeys.devSettings := Seq("play.server.http.idleTimeout" -> "220s")
