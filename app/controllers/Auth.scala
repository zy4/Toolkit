package controllers

import java.time.ZonedDateTime

import javax.inject.{ Inject, Singleton }
import actors.WebSocketActor.{ ChangeSessionID, LogOut }
import akka.actor.ActorRef
import de.proteinevolution.auth.UserSessions
import de.proteinevolution.auth.models.{ FormDefinitions, JSONTemplate }
import de.proteinevolution.models.database.users.{ User, UserToken }
import models.tools.ToolFactory
import de.proteinevolution.db.MongoStore
import de.proteinevolution.auth.models.MailTemplate._
import play.api.cache._
import play.api.i18n.I18nSupport
import play.api.mvc._
import play.api.libs.mailer._
import reactivemongo.bson._
import org.webjars.play.WebJarsUtil
import play.api.{ Environment, Logger }

import scala.concurrent.{ ExecutionContext, Future }

@Singleton
final class Auth @Inject()(
    webJarsUtil: WebJarsUtil,
    mongoStore: MongoStore,
    toolFactory: ToolFactory,
    userSessions: UserSessions,
    @NamedCache("wsActorCache") implicit val wsActorCache: SyncCacheApi,
    environment: Environment,
    assets: AssetsFinder,
    cc: ControllerComponents
)(implicit ec: ExecutionContext, mailerClient: MailerClient)
    extends AbstractController(cc)
    with I18nSupport
    with JSONTemplate
    with CommonController {

  private val logger = Logger(this.getClass)

  /**
   * Submission of the sign in form
   * Checks the Database for the user and logs him in if password matches
   *
   * @return
   */
  def signInSubmit: Action[AnyContent] =
    Action.async { implicit request =>
      userSessions.getUser.flatMap { unregisteredUser =>
        if (unregisteredUser.accountType < 0) {
          // Evaluate the Form
          FormDefinitions.signIn.bindFromRequest.fold(
            _ =>
              Future.successful {
                Ok(loginError())
            },
            // if no error, then insert the user to the collection
            signInFormUser => {
              val futureUser = mongoStore.findUser(
                BSONDocument(
                  "$or" -> List(BSONDocument(User.EMAIL -> signInFormUser.nameLogin),
                                BSONDocument(User.NAMELOGIN -> signInFormUser.nameLogin))
                )
              )
              futureUser.flatMap {
                case Some(databaseUser) =>
                  // Check the password
                  if (databaseUser.checkPassword(signInFormUser.password) && databaseUser.accountType > 0) {
                    // create a modifier document to change the last login date in the Database
                    val selector = BSONDocument(User.IDDB -> databaseUser.userID)
                    // Change the login time and give the new Session ID to the user.
                    // Additionally add the watched jobs to the users watchlist.
                    val modifier = userSessions.getUserModifier(databaseUser, forceSessionID = true)
                    // TODO this adds the non logged in user's jobs to the now logged in user's job list
                    //                            "$addToSet"        ->
                    //               BSONDocument(User.JOBS          ->
                    //               BSONDocument("$each"            -> unregisteredUser.jobs)))
                    // Finally add the edits to the collection
                    userSessions.modifyUserWithCache(selector, modifier).map {
                      case Some(loggedInUser) =>
                        logger.info(
                          "\n-[old user]-\n"
                          + unregisteredUser.toString
                          + "\n-[new user]-\n"
                          + loggedInUser.toString
                        )
                        // Remove the old, not logged in user
                        //removeUser(BSONDocument(User.IDDB -> unregisteredUser.userID))
                        userSessions.removeUserFromCache(unregisteredUser)

                        // Tell the job actors to copy all jobs connected to the old user to the new user
                        wsActorCache.get[List[ActorRef]](unregisteredUser.userID.stringify) match {
                          case Some(wsActors) =>
                            val actorList: List[ActorRef] = wsActors: List[ActorRef]
                            wsActorCache.set(loggedInUser.userID.stringify, actorList)
                            actorList.foreach(_ ! ChangeSessionID(loggedInUser.sessionID.get))
                            wsActorCache.remove(unregisteredUser.userID.stringify)
                          case None =>
                        }

                        // Everything is ok, let the user know that they are logged in now
                        Ok(loggedIn(loggedInUser))
                          .withSession(
                            userSessions.sessionCookie(request, loggedInUser.sessionID.get)
                          )
                      case None =>
                        Ok(loginIncorrect())
                    }
                  } else if (databaseUser.accountType < 1) {
                    // User needs to Verify first
                    Future.successful(Ok(mustVerify()))
                  } else {
                    // Wrong Password, show the error message
                    Future.successful(Ok(loginIncorrect()))
                  }
                case None =>
                  Future.successful {
                    Ok(loginIncorrect())
                  }
              }
            }
          )
        } else {
          Future.successful(Ok(alreadyLoggedIn()))
        }
      }
    }

  /**
   * Verifies a Token which was sent to the Users eMail address.
   * Token Types: 1 - eMail verification
   * 2 - password change verification
   * 3 - password reset verification
   * 4 -                             + reset
   *
   *
   * @param token
   * @return
   */
  def verification(nameLogin: String, token: String): Action[AnyContent] = Action.async { implicit request =>
    userSessions.getUser.flatMap { user: User =>
      // Grab the user from the database in case that the logged in user is not the user to verify
      // TODO check for both name or email
      mongoStore.findUser(BSONDocument(User.NAMELOGIN -> nameLogin)).flatMap {
        case Some(userToVerify) =>
          userToVerify.userToken match {
            case Some(userToken) =>
              if (userToken.token == token) {
                userToken.tokenType match {
                  case 1 => // Token for eMail verification
                    mongoStore
                      .modifyUser(
                        BSONDocument(User.IDDB -> userToVerify.userID),
                        BSONDocument(
                          "$set" ->
                          BSONDocument(User.ACCOUNTTYPE -> 1,
                                       User.DATEUPDATED -> BSONDateTime(ZonedDateTime.now.toInstant.toEpochMilli)),
                          BSONDocument(
                            "$unset" ->
                            BSONDocument(User.USERTOKEN -> "")
                          )
                        )
                      )
                      .map {
                        case Some(_) =>
                          Ok(
                            views.html.main(assets,
                                            webJarsUtil,
                                            toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                                            "Account verification was successful. Please log in.",
                                            "",
                                            environment)
                          )
                        case None => // Could not save the modified user to the DB
                          Ok(
                            views.html.main(
                              assets,
                              webJarsUtil,
                              toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                              "Verification was not successful due to a database error. Please try again later.",
                              "",
                              environment
                            )
                          )
                      }
                  case 2 => // Token for password change validation
                    userToken.passwordHash match {
                      case Some(newPassword) =>
                        mongoStore
                          .modifyUser(
                            BSONDocument(User.IDDB -> userToVerify.userID),
                            BSONDocument(
                              "$set" ->
                              BSONDocument(
                                User.PASSWORD    -> newPassword,
                                User.DATEUPDATED -> BSONDateTime(ZonedDateTime.now.toInstant.toEpochMilli)
                              ),
                              "$unset" ->
                              BSONDocument(User.SESSIONID -> "", User.CONNECTED -> "", User.USERTOKEN -> "")
                            )
                          )
                          .map {
                            case Some(modifiedUser) =>
                              userSessions.removeUserFromCache(user)
                              val eMail = PasswordChangedMail(modifiedUser)
                              eMail.send
                              // Force Log Out on all connected users.
                              (wsActorCache.get(modifiedUser.userID.stringify): Option[List[ActorRef]]) match {
                                case Some(webSocketActors) =>
                                  webSocketActors.foreach(_ ! LogOut)
                                case None =>
                              }
                              // User modified properly
                              Ok(
                                views.html.main(
                                  assets,
                                  webJarsUtil,
                                  toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                                  "Password change verification was successful. Please log in with Your new password.",
                                  "",
                                  environment
                                )
                              )
                            case None => // Could not save the modified user to the DB - failsave in case the DB is down
                              Ok(
                                views.html.main(
                                  assets,
                                  webJarsUtil,
                                  toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                                  "Verification was not successful due to a database error. Please try again later.",
                                  "",
                                  environment
                                )
                              )
                          }
                      case None =>
                        // This should not happen - Failsafe when the password hash got overwritten somehow
                        Future.successful(
                          Ok(
                            views.html
                              .main(
                                assets,
                                webJarsUtil,
                                toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                                "The Password you have entered was insufficient, please create a new one.",
                                "",
                                environment
                              )
                          )
                        )
                    }

                  case 3 =>
                    // Give a token to the current user to allow him to change the password in a different view (Password Recovery)
                    val newToken =
                      UserToken(tokenType = 4, token = userToken.token, userID = Some(userToVerify.userID))
                    val selector = BSONDocument(User.IDDB -> user.userID)
                    val modifier = BSONDocument(
                      "$set" -> BSONDocument(
                        User.DATEUPDATED -> BSONDateTime(ZonedDateTime.now.toInstant.toEpochMilli),
                        User.USERTOKEN   -> newToken
                      )
                    )
                    userSessions.modifyUserWithCache(selector, modifier).map {
                      case Some(_) =>
                        Ok(
                          views.html.main(assets,
                                          webJarsUtil,
                                          toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                                          "",
                                          "passwordReset",
                                          environment)
                        )
                      case None => // Could not save the modified user to the DB
                        Ok(
                          views.html.main(
                            assets,
                            webJarsUtil,
                            toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                            "Verification was not successful due to a database error. Please try again later.",
                            "",
                            environment
                          )
                        )
                    }
                  case _ =>
                    Future.successful(
                      Ok(
                        views.html.main(assets,
                                        webJarsUtil,
                                        toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                                        "There was an error finding your token.",
                                        "",
                                        environment)
                      )
                    )
                }

              } else {
                // No Token in DB
                Future.successful(
                  Ok(
                    views.html.main(assets,
                                    webJarsUtil,
                                    toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                                    "The token you used is not valid.",
                                    "",
                                    environment)
                  )
                )
              }
            case None =>
              Future.successful(
                Ok(
                  views.html.main(assets,
                                  webJarsUtil,
                                  toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                                  "There was an error finding your token.",
                                  "",
                                  environment)
                )
              )
          }
        case None =>
          Future.successful(
            Ok(
              views.html.main(assets,
                              webJarsUtil,
                              toolFactory.values.values.toSeq.sortBy(_.toolNameLong),
                              "There was an error finding your account.",
                              "",
                              environment)
            )
          )
      }
    }
  }

}
