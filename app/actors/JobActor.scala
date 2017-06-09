package actors
import javax.inject.{Inject, Named}
import play.api.Logger
import actors.JobActor._
import akka.actor.{Actor, ActorRef}
import akka.event.LoggingReceive
import models.Constants
import models.database.jobs._
import models.database.statistics.{JobEvent, JobEventLog}
import models.database.users.User
import models.mailing.JobFinishedMail
import models.search.JobDAO
import modules.tel.TEL
import modules.tel.runscripts._
import better.files._
import com.typesafe.config.ConfigFactory
import controllers.{FileException, UserSessions}
import models.sge.Qdel
import modules.{CommonModule, LocationProvider}
import modules.tel.env.Env
import modules.tel.execution.ExecutionContext.FileAlreadyExists
import modules.tel.execution.{ExecutionContext, RunningExecution, WrapperExecutionFactory}
import modules.tel.runscripts.Runscript.Evaluation
import org.joda.time.DateTime
import play.api.cache.{CacheApi, NamedCache}
import play.api.libs.mailer.MailerClient
import play.modules.reactivemongo.ReactiveMongoApi
import reactivemongo.bson.{BSONDateTime, BSONDocument, BSONObjectID}

import scala.concurrent.ExecutionContext.Implicits.global
import play.api.libs.json._

import scala.concurrent.Future
import scala.sys.process.Process
import scala.util.{Failure, Success}

object JobActor {

  case class PrepareJob(job: Job,
                        params: Map[String, String],
                        startJob: Boolean = false,
                        isInternalJob: Boolean = false)

  case class CheckJobHashes(jobID: String)

  // Messages the Job Actor to start a job
  case class StartJob(jobID: String)

  trait Factory {

    def apply: Actor
  }

  // Messages the jobActor accepts from outside
  case class PushJob(job: Job)

  // Message for the Websocket Actor to send a ClearJob Message
  case class ClearJob(jobID: String, deleted: Boolean = false)

  // Message for the Websocket Actor to send a ClearJob Message
  case class DeleteJob(jobID: String, userID: BSONObjectID)

  // User Actor starts watching
  case class AddToWatchlist(jobID: String, userID: BSONObjectID)

  // UserActor Stops Watching this Job
  case class RemoveFromWatchlist(jobID: String, userID: BSONObjectID)

  // JobActor is requested to Delete the job
  case class Delete(jobID: String, userID: BSONObjectID)

  // JobActor is requested to Delete the job
  case class DeleteFromDisk(job: Job)

  // JobActor is requested to mark all jobs that are older than
  // a given number of days
  case class MarkForDeletion(daysAgo: Int)

  // Job Controller receives a job state change from the SGE or from any other valid source
  case class JobStateChanged(jobID: String, jobState: JobState)

  // Job Controller receives push message to update the log
  case class UpdateLog(jobID: String)
}

class JobActor @Inject()(runscriptManager: RunscriptManager, // To get runscripts to be executed
                         env: Env, // To supply the runscripts with an environment
                         val reactiveMongoApi: ReactiveMongoApi,
                         implicit val mailerClient: MailerClient,
                         val jobDao: JobDAO,
                         qdel: Qdel,
                         wrapperExecutionFactory: WrapperExecutionFactory,
                         implicit val locationProvider: LocationProvider,
                         @Named("jobIDActor") jobIDActor: ActorRef,
                         @NamedCache("userCache") implicit val userCache: CacheApi,
                         @NamedCache("wsActorCache") implicit val wsActorCache: CacheApi)
    extends Actor
    with Constants
    with UserSessions
    with CommonModule {

  // Attributes asssocidated with a Job
  private var currentJobs: Map[String, Job]                           = Map.empty
  private var currentJobLogs: Map[String, JobEventLog]                = Map.empty
  private var currentExecutionContexts: Map[String, ExecutionContext] = Map.empty

  // Running executions
  private var runningExecutions: Map[String, RunningExecution] = Map.empty

  /**
    * Get the job from the current jobs. If it is not there, get it from the DB.
    *
    * @param jobID
    * @return
    */
  private def getCurrentJob(jobID: String): Future[Option[Job]] = {
    // Check if the job is still in the current jobs.
    this.currentJobs.get(jobID) match {
      case Some(job) => // Everything is fine. Return the job.
        Future.successful(Some(job))
      case None => // Job is not in the current jobs.. try to get it back.
        findJob(BSONDocument(Job.JOBID -> jobID)).map {
          case Some(job) =>
            // Get the job back into the current jobs
            this.currentJobs = this.currentJobs.updated(job.jobID, job)
            // TODO Check if the job is a running job and also if the cluster has done any changes with on the job.
            // Return the job
            Some(job)
          case None =>
            // Something must have went wrong, the job is not in the DB.
            None
        }
    }
  }

  /**
    * Gets the execution context for a given jobID, even if it has been removed.
    *
    * @param jobID
    * @return
    */
  private def getCurrentExecutionContext(jobID: String): Option[ExecutionContext] = {
    this.currentExecutionContexts.get(jobID) match {
      case Some(executionContext) => Some(executionContext)
      case None =>
        if ((jobPath / jobID).exists) {
          val executionContext = ExecutionContext(jobPath / jobID, reOpen = true)
          this.currentExecutionContexts = this.currentExecutionContexts.updated(jobID, executionContext)
          Some(executionContext)
        } else {
          None
        }
    }
  }

  /**
    * Return the validated parameters
    *
    * @param job
    * @param runscript
    * @param params
    * @return
    */
  private def validatedParameters(job: Job,
                                  runscript: Runscript,
                                  params: Map[String, String]): Seq[(String, (Evaluation, Option[Argument]))] = {
    // Representation of the current State of the job submission

    // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
    // TODO Check parameters for validity here!!!             TODO
    // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

    var validParameters: Seq[(String, (Evaluation, Option[Argument]))] =
      runscript.parameters.map(t => t._1 -> (t._2 -> Some(ValidArgument(new LiteralRepresentation(RString("false"))))))
    params.foreach { pv =>
      validParameters = supply(job.jobID, pv._1, pv._2, validParameters)
    }
    validParameters
  }

  /** Supplies a value for a particular Parameter. Returns params again if the parameter
    * is not present
    *
    * @param name
    * @param value
    * @param params
    */
  private def supply(jobID: String,
                     name: String,
                     value: String,
                     params: Seq[(String, (Runscript.Evaluation, Option[Argument]))])
    : Seq[(String, (Runscript.Evaluation, Option[Argument]))] = {
    params.map {
      case (paramName, (evaluation, _)) if paramName == name =>
        val x = Some(evaluation(RString(value), this.currentExecutionContexts(jobID)))
        (name, (evaluation, x))
      case q => q
    }
  }

  /**
    * JobActor removes the Job from its maps
    *
    * @param jobID
    */
  private def removeJob(jobID: String) = {
    val job = this.currentJobs(jobID)

    // If the job is in the current jobs remove it
    this.currentJobs = this.currentJobs.filter(_._1 == jobID)

    // Save Job Event Log to the collection and remove it from the map afterwards
    if (this.currentJobLogs.contains(jobID)) {
      addJobLog(this.currentJobLogs(jobID))
      this.currentJobLogs = this.currentJobLogs.filter(_._1 == jobID)
    }

    // If the job Appears in the running Execution, terminate it
    if (this.runningExecutions.contains(jobID)) {
      this.runningExecutions(jobID).terminate()
      this.runningExecutions = this.runningExecutions.filter(_._1 == jobID)
    }
    this.currentExecutionContexts = this.currentExecutionContexts.filter(_._1 == jobID)
  }

  /**
    * Trys to delete a job and inform all watching users about it
    *
    * @param job
    * @param userID
    */
  private def delete(job: Job, userID: BSONObjectID): Unit = {
    if (job.ownerID.contains(userID)) {
      Logger.info("Owner Requested job Deletion")

      // Message user clients to remove the job from their watchlist
      Logger.info("Informing Users of deletion.")
      val foundWatchers = job.watchList.flatMap(userID => wsActorCache.get(userID.stringify): Option[List[ActorRef]])
      foundWatchers.flatten.foreach(_ ! ClearJob(job.jobID))

      // Mark the job in mongoDB
      modifyJob(
        BSONDocument(Job.IDDB -> job.mainID),
        BSONDocument(
          "$set" ->
            BSONDocument(Job.DELETION -> JobDeletion(JobDeletionFlag.OwnerRequest, Some(DateTime.now()))),
          "$unset" ->
            BSONDocument(Job.WATCHLIST -> "")
        )
      ).foreach {
        case Some(deletedJob) =>
          Logger.info(s"Job Deletion from DB was successful:\n${deletedJob.toString()}")
        case None =>
          Logger.info("Job Deletion from DB failed.")
      }
    } else {
      // Just clear a job which is not owned by the user
      modifyJob(BSONDocument(Job.IDDB -> job.mainID), BSONDocument("$pull" -> BSONDocument(Job.WATCHLIST -> userID)))
    }
    // clear job from the user's watchlist
    modifyUserWithCache(BSONDocument(User.IDDB -> userID),
                        BSONDocument("$pull"   -> BSONDocument(User.JOBS -> job.jobID)))
  }

  /**
    * Updates Jobstate in Model, in database, and notifies user watchlist
    */
  private def updateJobState(job: Job): Future[Job] = {
    // Push the updated job into the current jobs
    this.currentJobs = this.currentJobs.updated(job.jobID, job)

    // Update job in the database and notify watcher upon completion
    modifyJob(BSONDocument(Job.IDDB -> job.mainID), BSONDocument("$set" -> BSONDocument(Job.STATUS -> job.status)))
      .map { modifiedJob =>
        val jobLog = this.currentJobLogs.get(job.jobID) match {
          case Some(jobEventLog) => jobEventLog.addJobStateEvent(job.status)
          case None =>
            JobEventLog(mainID = job.mainID,
                        toolName = job.tool,
                        events = List(JobEvent(job.status, Some(DateTime.now))))
        }
        this.currentJobLogs = this.currentJobLogs.updated(job.jobID, jobLog)
        val foundWatchers = job.watchList.flatMap(userID => wsActorCache.get(userID.stringify): Option[List[ActorRef]])
        foundWatchers.flatten.foreach(_ ! PushJob(job))
        job
      }
  }

  /**
    * Determines whether the parameter list is completely supplied
    *
    * @param params
    * @return
    */
  private def isComplete(params: Seq[(String, (Runscript.Evaluation, Option[Argument]))]): Boolean = {
    // If we have an argument for all parameters, we are done
    params.forall(item => item._2._2.isDefined)
  }

  /**
    * Sends an eMail to the owner of the job
    *
    * @param job
    * @return
    */
  private def sendJobUpdateMail(job: Job): Boolean = {
    if (job.emailUpdate && job.ownerID.isDefined) {
      findUser(BSONDocument(User.IDDB -> job.ownerID)).foreach {
        case Some(user) =>
          Logger.info("Sending EMail to owner: " + user.getUserData.eMail)
          val eMail = JobFinishedMail(user, job)
          eMail.send
        case None =>
      }
      true
    } else {
      false
    }
  }

  override def receive = LoggingReceive {

    case PrepareJob(job, params, startJob, isInternalJob) =>
      // jobid will also be available as parameter
      val extendedParams = params + ("jobid" -> job.jobID)

      // Add job to the current jobs
      this.currentJobs = this.currentJobs.updated(job.jobID, job)

      try {
        // Establish execution context for the new Job
        val executionContext = ExecutionContext(jobPath / job.jobID)
        this.currentExecutionContexts = this.currentExecutionContexts.updated(job.jobID, executionContext)

        // Create a log for this job
        this.currentJobLogs =
          this.currentJobLogs.updated(job.jobID,
            JobEventLog(mainID = job.mainID,
              toolName = job.tool,
              internalJob = isInternalJob,
              events = List(JobEvent(job.status, Some(DateTime.now)))))

        // Update the statistics
        increaseJobCount(job.tool) // TODO switch to better statistic handling

        // Get new runscript instance from the runscript manager
        val runscript: Runscript = runscriptManager(job.tool).withEnvironment(env)

        // Validate the Parameters right away
        var validParameters = this.validatedParameters(job, runscript, extendedParams)

        // adds the params of the disabled controls from formData, sets value of those to "false"
        validParameters.filterNot(pv => extendedParams.contains(pv._1)).foreach { pv =>
          extendedParams.+(pv._1 -> "false")
        }

        // Serialize the JobParameters to the JobDirectory
        // Store the extended Parameters in the working directory for faster reloading
        executionContext.writeParams(extendedParams)

        if (isComplete(validParameters)) {
          // When the user wants to force the job to start without job hash check, then this will jump right to prepared
          if (startJob) {
            val jobHash = JobHash.generateJobHash(job, params, env, jobDao)
            hashCollection.flatMap(_.insert(jobHash))
            self ! StartJob(job.jobID)
          } else {
            Logger.info("JobID " + job.jobID + " will now be hashed.")
            self ! CheckJobHashes(job.jobID)
          }
        } else {
          // TODO Implement Me. This specifies what the JobActor should do
          // TODO when not all parameters have been specified or when they are invalid
          Logger.error("[JobActor.PrepareJob] The job " + job.jobID + " has invalid or missing parameters.")
          self ! JobStateChanged(job.jobID, Error)
        }
      } catch {
        case FileAlreadyExists(msg) =>
          Logger.error(
            "[JobActor.PrepareJob] The directory for job " + job.jobID + " already exists\n" +
              "[JobActor.PrepareJob] Stopping job since it can not be retrieved by user.")
          self ! JobStateChanged(job.jobID, Error)
      }

    /**
      * Checks the jobHashDB for matches and generates one for the job if there are none.
      */
    case CheckJobHashes(jobID) =>
      Logger.info("JobID " + jobID + " will now be hashed.")
      this.getCurrentJob(jobID).foreach {
        case Some(job) =>
          this.getCurrentExecutionContext(jobID) match {
            case Some(executionContext) =>
              // Ensure that the jobID is not being hashed
              val params = executionContext.reloadParams
              val jobHash = JobHash.generateJobHash(job, params, env, jobDao)
              Logger.info("JobHash: " + jobHash.toString)
              // Match the hash
              jobDao.matchHash(jobHash).map { richSearchResponse =>
                Logger.info("Retrieved richSearchResponse")
                Logger.info("success: " + richSearchResponse.getHits.getHits.map(_.getId).mkString(", "))
                Logger.info("hits: " + richSearchResponse.totalHits)

                // Generate a list of hits and convert them into a list of future option jobs
                val mainIDs = richSearchResponse.getHits.getHits.toList.map { hit =>
                  BSONObjectID.parse(hit.getId).getOrElse(BSONObjectID.generate())
                }

                // Find the Jobs in the Database
                findJobs(BSONDocument(Job.IDDB -> BSONDocument("$in" -> mainIDs))).map { jobList =>
                  val foundMainIDs = jobList.map(_.mainID)

                  if (jobList.exists(_.status == Done)) {
                    Logger.info("JobID " + jobID + " is a duplicate.")
                    self ! JobStateChanged(job.jobID, Pending)
                  } else {
                    Logger.info("JobID " + jobID + " will now be started.")
                    hashCollection.flatMap(_.insert(jobHash))
                    self ! StartJob(job.jobID)
                  }
                }
              }
            case None =>
              Logger.error("[JobActor.CheckJobHashes] Could not recreate execution context for jobID " + jobID)
          }
        case None =>
          Logger.error("[JobActor.CheckJobHashes] Could not find the jobID " + jobID + " in the Cache or DB.")
      }

    /**
      * Checks everything and deletes the Job.
      * Includes: Removing the job from the cluster if the job is still in the current jobs
      * Creating the job deleted object and inserting it to the database
      * Removing the job from ES
      */
    case Delete(jobID, userID) =>
      Logger.info(s"Received Delete for $jobID")
      this.getCurrentJob(jobID).foreach {
        case Some(job) =>
          Logger.info("Removing Job from Elastic Search.")
          jobDao.deleteJob(job.mainID.stringify) // Remove job from elastic search
          Logger.info("Removing Job from current Jobs.")
          this.removeJob(jobID)
          Logger.info("Removing Job from DB")
          this.delete(job, userID)
          Logger.info("Deletion Complete.")
        case None =>
          Logger.info("No such jobID found in current jobs. Loading job from DB.")
          findJob(BSONDocument(Job.JOBID -> jobID)).map {
            case Some(job) =>
              Logger.info("Found Job in DB. Deleting.")
              Logger.info("Removing Job from Elastic Search.")
              jobDao.deleteJob(job.mainID.stringify) // Remove job from elastic search
              Logger.info("Removing Job from DB")
              this.delete(job, userID)
            case None =>
              Logger.error("[JobActor.Delete] No such jobID found in Database. Ignoring.")
          }
      }


    /**
      * Checks everything and deletes the Job from disk.
      * Includes: Removing the job from the cluster if the job is still in the current jobs
      * Creating the job deleted object and inserting it to the database
      * Removing the job from ES
      */
    case DeleteFromDisk(job) =>
      Logger.info(s"Received Delete from disk for ${job.jobID}")
      Logger.info("Removing Job from Elastic Search.")
      jobDao.deleteJob(job.mainID.stringify) // Remove job from elastic search
      Logger.info("Removing Job from current Jobs.")
      this.removeJob(job.jobID)
      Logger.info("Removing Job from mongo DB")
      removeJob(BSONDocument(Job.ID -> job.jobID))
      Logger.info("Deletion Complete.")
      val deletetionFile = s"$jobPath${job.jobID}/0/delete.sh".toFile
      if(!deletetionFile.exists) {
        throw FileException(s"File ${deletetionFile.name} does not exist.")
      }
      else if (!deletetionFile.isExecutable) {
        throw FileException(s"File ${deletetionFile.name} is not executable.")
      } else {
        Process(deletetionFile.pathAsString).run().exitValue() match {
          case 0 => Logger.info("Executing deletion.sh successful.")
          case _ => Logger.info("Executing deletion.sh failed.")
        }
      }


      /**
        * marks jobs in mongodb for deletion that are older than a given number of days
        * ('daysAgo') and informs all watching users
        * about it in behalf of the job maintenance routine
        *
        */
     case MarkForDeletion (daysAgo: Int) =>
       Logger.info("Removing Job from DB")
        modifyJob(
          BSONDocument(Job.DATECREATED -> BSONDocument( "$lt" -> BSONDateTime(new DateTime().minusDays(daysAgo).getMillis))),
          BSONDocument(
            "$set" ->
              BSONDocument(Job.DELETION -> JobDeletion(JobDeletionFlag.Automated, Some(DateTime.now()))),
            "$unset" ->
              BSONDocument(Job.WATCHLIST -> "")
          )
        ).foreach {
          case Some(deletedJob) =>
            Logger.info(s"Job mark for Deletion in DB was successful:\n${deletedJob.toString()}")
            Logger.info("Job maintenance routine Requested job Deletion")

            // Message user clients to remove the job from their watchlist
            Logger.info("Informing Users of deletion.")
            DeleteJob(deletedJob.jobID, deletedJob.ownerID.get)
            Logger.info("Removing Job from Elastic Search.")
            jobDao.deleteJob(deletedJob.mainID.stringify) // Remove job from elastic search


          case None =>
            Logger.info("Job mark Deletion in DB failed.")
        }

    /**
      * Starts the job
      */
    case StartJob(jobID) =>
      this.getCurrentJob(jobID).foreach {
        case Some(job) =>
          this.getCurrentExecutionContext(jobID) match {
            case Some(executionContext) =>
              Logger.info("[JobActor.StartJob] reached. starting job " + jobID)
              // set memory allocation on the cluster and let the clusterMonitor define the multiplier
              val queue = ConfigFactory.load().getString(s"Tools.${job.tool}.queue")
              val h_vmem = (ConfigFactory
                .load()
                .getString(s"Tools.${job.tool}.memory")
                .dropRight(1)
                .toInt * TEL.memFactor).toString + "G"
              val threads =
                math.ceil(ConfigFactory.load().getInt(s"Tools.${job.tool}.threads") * TEL.threadsFactor).toInt
              env.configure(s"MEMORY", h_vmem)
              env.configure(s"THREADS", threads.toString)
              env.configure(s"QUEUE", queue.toString)
              Logger.info(s"$jobID is running with $h_vmem h_vmem")
              Logger.info(s"$jobID is running with $threads threads")
              Logger.info(s"$jobID is queued: $queue")

              val clusterData = JobClusterData("", Some(h_vmem), Some(threads), Some(queue))

              modifyJob(BSONDocument(Job.IDDB -> job.mainID),
                        BSONDocument("$set" ->
                          BSONDocument(Job.CLUSTERDATA -> clusterData))).foreach {
                case Some(updatedJob) =>
                  // Get new runscript instance from the runscript manager
                  val runscript: Runscript = runscriptManager(job.tool).withEnvironment(env)
                  // Load the parameters from the serialized parameters file
                  val params = executionContext.reloadParams
                  // Validate the Parameters (again) to ensure that everything works
                  var validParameters = this.validatedParameters(job, runscript, params)

                  // adds the params of the disabled controls from formData, sets value of those to "false"
                  validParameters.filterNot(pv => params.contains(pv._1)).foreach { pv =>
                    params.+(pv._1 -> "false")
                  }

                  if (isComplete(validParameters)) {
                    val pendingExecution = wrapperExecutionFactory.getInstance(
                      runscript(validParameters.map(t => (t._1, t._2._2.get.asInstanceOf[ValidArgument]))))

                    if ( ! executionContext.blocked) {

                      executionContext.accept(pendingExecution)
                      Logger.info("[JobActor.StartJob] Running job now.")
                      this.runningExecutions = this.runningExecutions.updated(job.jobID, executionContext.executeNext.run())
                    }
                  } else {
                    // TODO Implement Me. This specifies what the JobActor should do if not all parameters have been specified
                    Logger.info("STAY")
                  }

                  self ! JobStateChanged(job.jobID, Prepared)
                case None =>
                  Logger.error("[JobActor.StartJob] Job could not be written to DB: " + jobID)
              }
            //env.remove(s"MEMORY")
            //env.remove(s"THREADS")
            case None =>
          }
        case None =>
          Logger.error("[JobActor.StartJob] Job not found in DB: " + jobID)
      }

    // User does no longer watch this Job (delete also from JobManager)
    case DeleteJob(jobID, userID) =>
      Logger.info("Removing Job from current Jobs.")
      this.removeJob(jobID)
      modifyJob(BSONDocument(Job.JOBID -> jobID), BSONDocument("$pull" -> BSONDocument(Job.WATCHLIST -> userID)))
        .foreach {
          case Some(updatedJob) =>
            this.currentJobs = this.currentJobs.updated(jobID, updatedJob)
          case None =>
        }
      val wsActors = wsActorCache.get(userID.stringify): Option[List[ActorRef]]
      wsActors.foreach(_.foreach(_ ! ClearJob(jobID, deleted = true)))

      modifyUserWithCache(BSONDocument(User.IDDB -> userID), BSONDocument("$pull" -> BSONDocument(User.JOBS -> jobID)))

    // User Starts watching job
    case AddToWatchlist(jobID, userID) =>
      modifyJob(BSONDocument(Job.JOBID -> jobID), BSONDocument("$addToSet" -> BSONDocument(Job.WATCHLIST -> userID)))
        .map {
          case Some(updatedJob) =>
            modifyUserWithCache(BSONDocument(User.IDDB   -> userID),
                                BSONDocument("$addToSet" -> BSONDocument(User.JOBS -> jobID))).foreach { _ =>
              this.currentJobs = this.currentJobs.updated(jobID, updatedJob)
              val wsActors = wsActorCache.get(userID.stringify): Option[List[ActorRef]]
              wsActors.foreach(_.foreach(_ ! PushJob(updatedJob)))
            }
          case None =>
        }

    // User does no longer watch this Job (stays in JobManager)
    case RemoveFromWatchlist(jobID, userID) =>
      modifyJob(BSONDocument(Job.JOBID -> jobID), BSONDocument("$pull" -> BSONDocument(Job.WATCHLIST -> userID)))
        .foreach {
          case Some(updatedJob) =>
            modifyUserWithCache(BSONDocument(User.IDDB -> userID),
                                BSONDocument("$pull"   -> BSONDocument(User.JOBS -> jobID))).foreach { _ =>
              this.currentJobs = this.currentJobs.updated(jobID, updatedJob)
              val wsActors = wsActorCache.get(userID.stringify): Option[List[ActorRef]]
              wsActors.foreach(_.foreach(_ ! ClearJob(jobID)))
            }
          case None =>
        }

    // Message from outside that the jobState has changed
    case JobStateChanged(jobID: String, jobState: JobState) =>
      this.getCurrentJob(jobID).foreach {
        case Some(oldJob) =>
          // Update the job object
          val job = oldJob.copy(status = jobState)
          // Give a update message to all
          Logger.info("Job State has changed to " + job.status.toString + " for the Job with the JobID " + job.jobID)

          // Dependent on the state, we have to do different things
          job.status match {
            case Done =>
              // Job is no longer running
              Logger.info("Removing execution context")
              this.runningExecutions = this.runningExecutions.-(job.jobID)
              Logger.info("DONE Removing execution context")
              // Tell the user that their job finished via eMail
              sendJobUpdateMail(job)
              val result = (jobPath / job.jobID / "results").list
                .withFilter(_.hasExtension)
                .withFilter(_.extension.get == ".json")
                .map { file =>
                  (file.nameWithoutExtension,
                   reactivemongo.play.json.BSONFormats.toBSON(Json.parse(file.contentAsString)).get)
                }
                .toTraversable
              if(result.nonEmpty){
                // Put the result files into the database, JobActor has to wait until this process has finished
                val x = result2Job(job.jobID, BSONDocument(result)) onComplete {
                  case Success(doc) =>
                    // Now we can update the JobState and remove it, once the update has completed
                    this.updateJobState(job).map { job =>
                      this.removeJob(job.jobID)
                      Logger.info("Job has been removed from JobActor")
                    }
                  case Failure(t) => println("An error has occured: " + t.getMessage)
                }
              }else {
                // Now we can update the JobState and remove it, once the update has completed
                this.updateJobState(job).map { job =>
                  this.removeJob(job.jobID)
                  Logger.info("Job has been removed from JobActor")
                }
              }

            // Currently no further error handling
            case Error =>
              this.updateJobState(job).map { job =>
                this.removeJob(job.jobID)

                // Tell the user that their job failed via eMail
                sendJobUpdateMail(job)

                // Update the statistics for the failed job TODO - swap to better statistic handling
                increaseJobCount(job.tool, failed = true)
              }

            case _ =>
              this.updateJobState(job)
          }
        case None =>
          Logger.info("Job not found: " + jobID)
      }

    case UpdateLog(jobID: String) =>
      currentJobs.get(jobID) match {
        case Some(job) =>
          val foundWatchers =
            job.watchList.flatMap(userID => wsActorCache.get(userID.stringify): Option[List[ActorRef]])
          foundWatchers.flatten.foreach(_ ! PushJob(job))
        case None =>
      }
  }
}
