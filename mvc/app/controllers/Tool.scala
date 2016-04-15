package controllers

import javax.inject.{Named, Singleton, Inject}

import akka.util.Timeout
import models.tools._
import models.sessions.Session
import scala.concurrent.duration._
import play.api.i18n.{I18nSupport, MessagesApi}
import play.api.mvc.{Action, Controller}
import actors.MasterConnection

/**
  * Created by lukas on 1/27/16.
  */
object Tool {

  val tools : List[ToolModel] = List(Hmmer3, Tcoffee, Alnviz, Psiblast) // list of all added tools


  /** getToolModel
    * Returns the tool object for a tool's name, null when there is no such tool.
    *
    * @param toolName tool Name
    * @return
    */
  def getToolModel(toolName: String): ToolModel = {
    for (tool <- tools) {
      if ( tool.toolNameShort        == toolName
        || tool.toolNameLong         == toolName
        || tool.toolNameAbbreviation == toolName)  return tool
    }
    null
  }
}

@Singleton
class Tool @Inject()(val messagesApi: MessagesApi,
                     masterConnection: MasterConnection) extends Controller with I18nSupport {

  implicit val timeout = Timeout(5.seconds)

  import models.distributed.FrontendMasterProtocol._

  def submit(toolname: String, start : Boolean, newJob : Boolean) = Action { implicit request =>

    val sessionID = Session.requestSessionID(request) // Grab the Session ID

    // Fetch the job ID from the submission, might be the empty string
    val jobID = request.body.asFormUrlEncoded.get("jobid").head


    // Determine whether the toolname was fine
    val tool = Tool.getToolModel(toolname)

    // Check if the tool name was ok.
    if (tool == null)  NotFound
    else {

      // TODO replace with reflection to avoid the need to mention each tool explicitly here
      val form = tool.toolNameShort match {
        case "alnviz" => Alnviz.inputForm
        case "tcoffee" => Tcoffee.inputForm
        case "hmmer3" => Hmmer3.inputForm
        case "psiblast" => Psiblast.inputForm
        case "reformat" => Reformat.inputForm
      }
      val boundForm = form.bindFromRequest

      boundForm.fold(
        formWithErrors => {

          BadRequest("This was an error")
        },
        _ => {

          if(start) {

            masterConnection.masterProxy ! PrepareAndStart(sessionID, jobID, toolname, boundForm.data, newJob)
          } else {

            masterConnection.masterProxy ! Prepare(sessionID, jobID, toolname, boundForm.data, newJob)
          }
        }
      )
      Ok
    }
  }
}

/*
case class Prepare(sessionID : String,
                     jobID : String,
                     toolname : String,
                     params : Map[String, String],
                     newJob : Boolean) extends UserRequestWithJob(sessionID, jobID)

 */


/*
// Determine the submitted JobID
    val job_id =  request.body.asFormUrlEncoded.get("jobid").head match {

      // If the user has not provided any, this will be None
      case m if m.isEmpty => None

      // If the user has provided one or the Job is an already exisiting one, then there is a Job ID
      case m => Some(m)
    }

    Logger.info("Submission for JobID " + job_id.toString + " received")

    (userManager ? GetUserActor(session_id)).mapTo[ActorRef].map { userActor =>

      val tool = Tool.getToolModel(toolname)
      // Check if the tool name was ok.
      if (tool == null)
        NotFound
      else {
        //val form = tool.inputForm
        // TODO replace with reflection
        val form = tool.toolNameShort match {
          case "alnviz" => Alnviz.inputForm
          case "tcoffee" => Tcoffee.inputForm
          case "hmmer3" => Hmmer3.inputForm
          case "psiblast" => Psiblast.inputForm
          case "reformat" => Reformat.inputForm
        }
        val boundForm = form.bindFromRequest

        boundForm.fold(
          formWithErrors => {
            BadRequest("This was an error")
          },
          _ => {

            userActor ! PrepWD(toolname, boundForm.data, startImmediate, job_id, newSubmission)
          }
        )
        Ok
      }
    }


 */
