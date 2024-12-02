import play.sbt.PlayImport._
import sbt._

object Dependencies {

  val akkaVersion = "2.6.20"
  val circeV      = "0.14.2"

  lazy val commonDeps: Seq[ModuleID] = Seq(
    ws,
    filters,
    ehcache,
    guice,
    "com.typesafe.akka"    %% "akka-actor"                 % akkaVersion,
    "com.typesafe.akka"    %% "akka-cluster"               % akkaVersion,
    "com.typesafe.akka"    %% "akka-cluster-tools"         % akkaVersion,
    "com.typesafe.akka"    %% "akka-cluster-metrics"       % akkaVersion,
    "com.typesafe.akka"    %% "akka-slf4j"                 % akkaVersion,
    "com.typesafe.akka"    %% "akka-stream"                % akkaVersion,
    "com.typesafe.akka"    %% "akka-serialization-jackson" % akkaVersion,
    "com.typesafe.akka"    %% "akka-actor-typed"           % akkaVersion,
    "com.typesafe.play"    %% "play-mailer"                % "8.0.1",
    "com.typesafe.play"    %% "play-mailer-guice"          % "8.0.1",
    "com.github.pathikrit" %% "better-files"               % "3.9.1",
    "org.mindrot"           % "jbcrypt"                    % "0.4",
    ("org.reactivemongo"   %% "play2-reactivemongo"        % "1.0.10-play28")
      .exclude("com.typesafe.akka", "*")  // provided
      .exclude("com.typesafe.play", "*"), // provided
    ("org.reactivemongo" %% "reactivemongo-akkastream" % "1.0.10")
      .exclude("com.typesafe.akka", "*")  // provided
      .exclude("com.typesafe.play", "*"), // provided
    "org.typelevel" %% "cats-core"            % "2.8.0",
    "org.typelevel" %% "cats-effect"          % "3.3.14",
    "co.fs2"        %% "fs2-core"             % "3.3.0",
    "co.fs2"        %% "fs2-io"               % "3.3.0",
    "com.chuusai"   %% "shapeless"            % "2.3.10",
    "com.dripower"  %% "play-circe"           % "2814.2",
    "io.circe"      %% "circe-generic"        % circeV,
    "io.circe"      %% "circe-generic-extras" % circeV,
    ("io.circe"     %% "circe-bson" % "0.5.0").exclude("org.reactivemongo", "*"), // provided by play2-reactivemongo
    "org.typelevel" %% "simulacrum" % "1.0.0",
    "com.github.dalet-oss" % "mongobee" % "1.0.4",
    "com.maxmind.geoip2"   % "geoip2"   % "2.13.0"
  )

  lazy val testDeps: Seq[ModuleID] = Seq(
    "com.typesafe.akka"        %% "akka-testkit"        % akkaVersion % Test,
    "com.typesafe.akka"        %% "akka-stream-testkit" % akkaVersion % Test,
    "org.scalatestplus.play"   %% "scalatestplus-play"  % "5.1.0"     % Test,
    "org.mockito"               % "mockito-core"        % "4.8.1"     % Test,
    "com.softwaremill.macwire" %% "macros"              % "2.5.9"     % Test,
    "org.awaitility"            % "awaitility"          % "4.2.0"     % Test
  )

}
