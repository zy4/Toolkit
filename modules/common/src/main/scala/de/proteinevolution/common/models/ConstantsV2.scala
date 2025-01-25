/*
 * Copyright 2018 Dept. of Protein Evolution, Max Planck Institute for Biology
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

package de.proteinevolution.common.models

import java.io.File

import javax.inject.{ Inject, Singleton }
import play.api.Configuration

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.matching.Regex

@Singleton
class ConstantsV2 @Inject() (config: Configuration) {

  /** Number of jobActors */
  final val nJobActors: Int = 100

  /** Key for the Modeller tool */
  final val modellerKey: String = config.get[String]("modeller_key")

  // File Settings
  /** File Seperator */
  final val SEPARATOR: String = File.separator

  /** Path to the jobs folder */
  final val jobPath: String = s"${config.get[String]("job_path")}$SEPARATOR"

  /** results folder name */
  final val resultFolderName: String = "results"

  /** name of the parameter file */
  final val serializedParam: String = "sparam"

  final val formMultiValueSeparator: String = " "

  // Job limitation settings
  final val maxJobNum: Int        = 200  // max number of jobs that can be submitted from one ip within maxJobsWithin
  final val maxJobsWithin: Int    = 60   // time in minutes within the max number of jobs is applied
  final val maxJobNumDay: Int     = 2000 // max number of jobs that can be submitted from one ip within maxJobsWithinDay
  final val maxJobsWithinDay: Int = 1    // time in days within the max number of jobs is applied for a day

  // Job deletion settings
  /** Sweeps at this time after server start */
  final val jobDeletionDelay: FiniteDuration = 10.minutes

  /** Sweeps in this interval */
  final val jobDeletionInterval: FiniteDuration = 3.hours

  /** all jobs of registered users that are older than the given number are permanently deleted everywhere */
  final val jobDeletionRegistered: Int = 90

  /** all jobs of non registered users that are older than the given number are permanently deleted everywhere */
  final val jobDeletion: Int = 24

  /** time in days in which the job should not have been viewed */
  final val jobDeletionLastViewed: Int = 7

  // User Deletion settings
  /** Sweeps at this time after server start */
  final val userDeletionDelay: FiniteDuration = 70.minutes

  /** Sweeps in this interval */
  final val userDeletionInterval: FiniteDuration = 3.hours

  /** Sending an eMail to the user this many days before the deletion */
  final val userDeletionWarning: Int = 14 // days
  /** Deletes regular accounts after this time frame */
  final val userDeleting: Int = 1 // months
  /** Deletes users awaiting registration after this time frame */
  final val userDeletingRegisterEmail: Int = 3 // days
  /** Deletes registered accounts after this time frame */
  final val userDeletingRegistered: Int = 24 // months

  // qstat long polling in order to detect zombie jobs, not used for cluster load calc anymore
  /** Interval of the qstat requests */
  final val pollingInterval: FiniteDuration = 5.minutes

  /** The marker for 100% load capacity */
  final val loadPercentageMarker: Int = 32 // Jobs

  // jobID pattern settings
  /** allowed elements in the jobID */
  final val jobIDCharacters: String = "[0-9a-zA-Z_]"

  /** versioning character */
  final val jobIDVersioningCharacter: String = "_"

  /** The regular jobID pattern to match against */
  final val jobIDNoVersionPattern: Regex = s"($jobIDCharacters{3,96})".r

  /** The additional pattern for versioning */
  final val jobVersionPattern: Regex = s"(?:$jobIDVersioningCharacter([0-9]{1,3}))".r

  /** The combined pattern */
  final val jobIDPattern: Regex = (jobIDNoVersionPattern.regex + jobVersionPattern.regex).r

  /** The combined pattern with the version as an option */
  final val jobIDVersionOptionPattern: Regex = (jobIDNoVersionPattern.regex + jobVersionPattern.regex + "?").r

}
