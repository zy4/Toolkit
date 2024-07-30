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

package de.proteinevolution.jobs.models

import java.time.ZonedDateTime

import de.proteinevolution.common.models.database.jobs.JobState
import de.proteinevolution.common.models.database.jobs.JobState._
import de.proteinevolution.common.models.util.ZonedDateTimeHelper
import de.proteinevolution.tools.ToolConfig
import de.proteinevolution.user.User
import io.circe.JsonObject
import io.circe.syntax._
import play.api.Configuration
import reactivemongo.api.bson._

import cats.effect.unsafe.implicits.global

case class Job(
    jobID: String,
    ownerID: String,
    parentID: Option[String] = None,
    hash: Option[String] = None,
    isPublic: Boolean = false,
    status: JobState = Submitted,
    emailUpdate: Boolean = false,
    tool: String,
    watchList: List[String] = List.empty,
    clusterData: Option[JobClusterData] = None,
    dateCreated: ZonedDateTime = ZonedDateTime.now,
    dateUpdated: ZonedDateTime = ZonedDateTime.now,
    dateViewed: ZonedDateTime = ZonedDateTime.now,
    dateDeletionOn: ZonedDateTime,
    IPHash: Option[String]
) {

  def jsonPrepare(
      toolConfig: ToolConfig,
      user: User,
      paramValues: Option[Map[String, String]] = None
  )(implicit config: Configuration): JsonObject = {
    val toolObj = toolConfig.values.unsafeRunSync().apply(tool)
    JsonObject(
      "jobID"        -> jobID.asJson,
      "parentID"     -> parentID.asJson,
      "status"       -> status.asJson,
      "foreign"      -> !ownerID.equals(user.userID).asJson,
      "watched"      -> watchList.contains(user.userID).asJson,
      "isPublic"     -> isPublic.asJson,
      "code"         -> toolObj.code.asJson,
      "tool"         -> tool.asJson,
      "toolnameLong" -> toolObj.toolNameLong.asJson,
      "paramValues"  -> paramValues.asJson,
      "views"        -> toolObj.resultViews.asJson,
      "dateCreated"  -> dateCreated.toInstant.toEpochMilli.asJson,
      "dateUpdated"  -> dateUpdated.toInstant.toEpochMilli.asJson,
      "dateViewed"   -> dateViewed.toInstant.toEpochMilli.asJson
    )
  }

  override def toString: String = {
    s"""--[Job Object]--
       |jobID: ${this.jobID}
       |parentID: ${this.parentID}
       |tool: ${this.tool}
       |state: ${this.status}
       |ownerID: ${this.ownerID}
       |created on: ${this.dateCreated.toString}
       |--[Job Object end]--
     """.stripMargin
  }

  def isFinished: Boolean = status == Done || status == Error

}

object Job {

  final val ID               = "id"
  final val PARENT_ID        = "parentID"
  final val HASH             = "hash"
  final val OWNER_ID         = "ownerID"
  final val IS_PUBLIC        = "isPublic"
  final val STATUS           = "status"
  final val EMAIL_UPDATE     = "emailUpdate"
  final val TOOL             = "tool"
  final val CODE             = "code"
  final val LABEL            = "label"
  final val WATCH_LIST       = "watchList"
  final val CLUSTER_DATA     = "clusterData"
  final val SGE_ID           = s"$CLUSTER_DATA.${JobClusterData.SGE_ID}"
  final val DATE_CREATED     = "dateCreated"
  final val DATE_UPDATED     = "dateUpdated"
  final val DATE_VIEWED      = "dateViewed"
  final val DATE_DELETION_ON = "dateDeletionOn"
  final val TOOLNAME_LONG    = "toolnameLong"
  final val IP_HASH          = "ipHash"

  // TODO Bson macros handler
  implicit def reader: BSONDocumentReader[Job] =
    BSONDocumentReader[Job] { bson =>
      Job(
        jobID = bson.getAsTry[String](ID).getOrElse("Error loading Job Name"),
        parentID = bson.getAsOpt[String](PARENT_ID),
        hash = bson.getAsOpt[String](HASH),
        ownerID = bson.getAsTry[String](OWNER_ID).getOrElse("Error loading Job Owner"),
        isPublic = bson.getAsTry[Boolean](IS_PUBLIC).getOrElse(false),
        status = bson.getAsTry[JobState](STATUS).getOrElse(Error),
        emailUpdate = bson.getAsTry[Boolean](EMAIL_UPDATE).getOrElse(false),
        tool = bson.getAsTry[String](TOOL).getOrElse(""),
        watchList = bson.getAsTry[List[String]](WATCH_LIST).getOrElse(List.empty),
        clusterData = bson.getAsOpt[JobClusterData](CLUSTER_DATA),
        dateCreated =
          bson.getAsTry[BSONDateTime](DATE_CREATED).map(ZonedDateTimeHelper.getZDT).getOrElse(ZonedDateTime.now),
        dateUpdated =
          bson.getAsTry[BSONDateTime](DATE_UPDATED).map(ZonedDateTimeHelper.getZDT).getOrElse(ZonedDateTime.now),
        dateViewed =
          bson.getAsTry[BSONDateTime](DATE_VIEWED).map(ZonedDateTimeHelper.getZDT).getOrElse(ZonedDateTime.now),
        dateDeletionOn = bson
          .getAsTry[BSONDateTime](DATE_DELETION_ON)
          .map(ZonedDateTimeHelper.getZDT)
          .getOrElse(ZonedDateTime.now.plusDays(5)),
        IPHash = bson.getAsOpt[String](IP_HASH)
      )
    }

  implicit def writer: BSONDocumentWriter[Job] =
    BSONDocumentWriter[Job] { job =>
      BSONDocument(
        ID               -> job.jobID,
        PARENT_ID        -> job.parentID,
        HASH             -> job.hash,
        OWNER_ID         -> job.ownerID,
        IS_PUBLIC        -> job.isPublic,
        STATUS           -> job.status,
        EMAIL_UPDATE     -> job.emailUpdate,
        TOOL             -> job.tool,
        WATCH_LIST       -> job.watchList,
        CLUSTER_DATA     -> job.clusterData,
        DATE_CREATED     -> BSONDateTime(job.dateCreated.toInstant.toEpochMilli),
        DATE_UPDATED     -> BSONDateTime(job.dateUpdated.toInstant.toEpochMilli),
        DATE_VIEWED      -> BSONDateTime(job.dateViewed.toInstant.toEpochMilli),
        DATE_DELETION_ON -> BSONDateTime(job.dateDeletionOn.toInstant.toEpochMilli),
        IP_HASH          -> job.IPHash
      )
    }
}
