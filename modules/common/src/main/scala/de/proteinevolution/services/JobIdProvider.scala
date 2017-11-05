package de.proteinevolution.services

import javax.inject.{ Inject, Singleton }
import de.proteinevolution.db.MongoStore
import scala.collection.mutable.ListBuffer
import scala.concurrent.{ Await, ExecutionContext, Future }
import scala.util.Random

@Singleton
class JobIdProvider @Inject()(mongoStore: MongoStore)(
    implicit ec: ExecutionContext
) {

  // don't use ids which are provided but not in the database yet
  private var usedIds = new ListBuffer[String]()

  private def isValid(id: String): Future[Boolean] = {
    mongoStore.selectJob(id).map(_.isEmpty)
  }

  def provide: String = {
    val id = Iterator
      .continually[String](Random.nextInt(9999999).toString.padTo(7, '0'))
      .filter(x => Await.result(isValid(x), scala.concurrent.duration.Duration.Inf))
      .filterNot(usedIds.contains)
      .next()
    usedIds += id
    id
  }

  def trash(id: String): Unit = {
    this.usedIds -= id
  }

  // TODO add mongodb constraint that jobIds is a unique field

}
