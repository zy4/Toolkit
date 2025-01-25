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

package de.proteinevolution.jobs.controllers
import de.proteinevolution.auth.services.UserSessionService
import de.proteinevolution.auth.util.UserAction
import de.proteinevolution.base.controllers.ToolkitController
import de.proteinevolution.common.models.ConstantsV2
import de.proteinevolution.jobs.actors.JobActor.{ CheckIPHash, Delete, SetJobPublic }
import de.proteinevolution.jobs.dao.JobDao
import de.proteinevolution.jobs.services._
import de.proteinevolution.tools.ToolConfig
import io.circe.syntax._
import io.circe.{ Json, JsonObject }
import javax.inject.{ Inject, Singleton }
import play.api.Logging
import play.api.mvc.{ Action, AnyContent, ControllerComponents }

import scala.concurrent.ExecutionContext

@Singleton
class SubmissionController @Inject() (
    jobActorAccess: JobActorAccess,
    userSessions: UserSessionService,
    jobDispatcher: JobDispatcher,
    constants: ConstantsV2,
    cc: ControllerComponents,
    jobDao: JobDao,
    toolConfig: ToolConfig,
    jobResubmitService: JobResubmitService,
    jobIdProvider: JobIdProvider,
    jobFrontendToolsService: JobFrontendToolsService,
    userAction: UserAction
)(implicit ec: ExecutionContext)
    extends ToolkitController(cc)
    with Logging {

  def startJob(jobID: String): Action[AnyContent] = userAction { implicit request =>
    jobActorAccess.sendToJobActor(jobID, CheckIPHash(jobID))
    Ok(JsonObject("message" -> Json.fromString("Starting Job...")).asJson)
  }

  def changeJob(jobID: String): Action[Json] = userAction(circe.json).async { implicit request =>
    jobDao.findJob(jobID).map {
      case Some(job) =>
        if (!job.ownerID.equals(request.user.userID)) {
          Unauthorized
        } else {
          request.body.asObject match {
            case None => BadRequest
            case Some(obj) =>
              if (obj.contains("isPublic")) {
                jobActorAccess.sendToJobActor(
                  jobID,
                  SetJobPublic(jobID, obj("isPublic").get.asBoolean.getOrElse(false))
                )
              }
              Ok
          }
        }
      case None => // job does not exist
        NotFound
    }
  }

  def frontend(toolName: String): Action[AnyContent] = Action.async { implicit request =>
    if (toolConfig.isTool(toolName)) {
      // Add Frontend Job to Database
      jobFrontendToolsService.logFrontendJob(toolName).map(_ => NoContent)
    } else {
      fuccess(BadRequest)
    }
  }

  def delete(jobID: String): Action[AnyContent] = userAction { implicit request =>
    logger.info("Delete Action in JobController reached")
    jobActorAccess.sendToJobActor(jobID, Delete(jobID, request.user.userID))
    NoContent
  }

  def submitJob(toolName: String): Action[Json] = userAction(circe.json).async { implicit request =>
    request.body.asObject match {
      case None => fuccess(BadRequest)
      case Some(obj) =>
        val parts: Iterable[(String, String)] = for {
          (key, json) <- obj.toIterable
          str = json.fold[String](
            "",
            bool => bool.toString,
            num => num.toString,
            identity,
            vec => vec.toString,
            obj => obj.toString
          )
        } yield (key, str)
        jobDispatcher
          .submitJob(
            toolName,
            parts.toMap,
            request.user
          )
          .value
          .map {
            case Right(job) =>
              Ok(
                JsonObject(
                  "successful" -> Json.fromBoolean(true),
                  "code"       -> Json.fromInt(0),
                  "message"    -> Json.fromString("Submission successful."),
                  "jobID"      -> Json.fromString(job.jobID)
                ).asJson
              )
            case Left(error) => BadRequest(errors(error.msg))
          }
    }
  }

  def checkJobID(newJobID: String): Action[AnyContent] = Action.async { implicit request =>
    jobResubmitService.checkJobID(newJobID).map(r => Ok(r.asJson))
  }

}
