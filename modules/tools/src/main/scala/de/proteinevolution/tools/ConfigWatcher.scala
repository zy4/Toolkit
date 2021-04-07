/*
 * Copyright 2018 Dept. Protein Evolution, Max Planck Institute for Developmental Biology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package de.proteinevolution.tools

import de.proteinevolution.tel.param.ParamCollector
import play.api.{Configuration, Logging}
import better.files._
import cats.effect.{IO, Resource, Ref}
import fs2.Stream
import fs2.io.file.Files
import fs2.io.Watcher

import cats.effect.unsafe.implicits.global

import java.nio.file.Path
import javax.inject.{Inject, Singleton}

@Singleton
final private[tools] class ConfigWatcher @Inject()(
    pc: ParamCollector,
    config: Configuration,
    toolConfig: ToolConfig
) extends Logging {

  private[this] final val REFRESH_FILE = "tel.params_refresh"

  // start the file watcher
  Stream
    .eval(for {
      refreshFile <- IO.fromOption(config.get[Option[String]](REFRESH_FILE))(
        new IllegalArgumentException(
          s"configured file $REFRESH_FILE is missing"))
      path = refreshFile.toFile.path
      _ <- IO(logger.info(s"using $path as trigger for param reload"))
      ref <- toolConfig.ref
    } yield (path, ref))
    .evalMap {
      case (p, r) =>
        watch(Resource.eval(IO(p)), r).repeat.compile.drain
    }
    .compile
    .drain
    .unsafeRunAsync(_ => ())

  // fs2 file watcher
  private[this] def watch(
      source: Resource[IO, Path],
      r: Ref[IO, Map[String, Tool]]
  ): Stream[IO, Unit] =
    Stream
      .resource(source)
      .flatMap { f =>
        Files[IO]
          .watch(f)
          .map {
            case Watcher.Event.Modified(_, _) | Watcher.Event.Created(_, _) =>
              logger.info(
                s"file $REFRESH_FILE changed, reloading parameters ...")
              r.modify(_ => (toolConfig.readFromFile(), ())) // update since the config has to be re-parsed
              pc.reloadValues() // reload params
            case Watcher.Event.Deleted(_, _) =>
              logger.warn(s"file $REFRESH_FILE was deleted")
            case Watcher.Event.Overflow(_) =>
              logger.warn(s"file $REFRESH_FILE overflow")
            case Watcher.Event.NonStandard(_, _) =>
              logger.warn(s"file $REFRESH_FILE changed unexpectedly")
          }
      }

}
