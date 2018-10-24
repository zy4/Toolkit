package de.proteinevolution.models.database.jobs

import reactivemongo.bson.{ BSONInteger, BSONReader, BSONWriter }

sealed trait JobState {
  def toInt: Int
}

object JobState {

  import io.circe.{ Decoder, HCursor }

  case object Prepared extends JobState {
    override def toInt = 1
  }

  case object Queued extends JobState {
    override def toInt = 2
  }

  case object Running extends JobState {
    override def toInt = 3
  }

  case object Error extends JobState {
    override def toInt = 4
  }

  case object Done extends JobState {
    override def toInt = 5
  }

  case object Submitted extends JobState {
    override def toInt = 6
  }

  case object Pending extends JobState {
    override def toInt = 7
  }

  case object LimitReached extends JobState {
    override def toInt = 8
  }

  case object Deleted extends JobState {
    override def toInt = 9
  }

  import de.proteinevolution.base.helpers.ToolkitTypes._
  import io.circe.Encoder
  import shapeless._

  implicit val jobStateDecoder: Decoder[JobState] = (c: HCursor) =>
    for {
      number <- c.downField("status").as[Int]
    } yield {
      implicitly[AllSingletons[
        JobState,
        Prepared.type :+:
        Queued.type :+:
        Running.type :+:
        Done.type :+:
        Error.type :+:
        Submitted.type :+:
        Pending.type :+:
        LimitReached.type
        :+: CNil
      ]].values.find(_.toInt == number).getOrElse(throw new Exception)
  }

  implicit val jobStateEncoder: Encoder[JobState] = Encoder[Int].contramap(_.toInt)

  implicit object JobStateReader extends BSONReader[BSONInteger, JobState] {
    def read(state: BSONInteger): JobState with Product with Serializable = {
      state match {
        case BSONInteger(1) => Prepared
        case BSONInteger(2) => Queued
        case BSONInteger(3) => Running
        case BSONInteger(4) => Error
        case BSONInteger(5) => Done
        case BSONInteger(6) => Submitted
        case BSONInteger(7) => Pending
        case BSONInteger(8) => LimitReached
        case BSONInteger(9) => Deleted
        case _              => Error
      }
    }
  }

  implicit object JobStateWriter extends BSONWriter[JobState, BSONInteger] {
    def write(state: JobState): BSONInteger = {
      state match {
        case Prepared     => BSONInteger(1)
        case Queued       => BSONInteger(2)
        case Running      => BSONInteger(3)
        case Error        => BSONInteger(4)
        case Done         => BSONInteger(5)
        case Submitted    => BSONInteger(6)
        case Pending      => BSONInteger(7)
        case LimitReached => BSONInteger(8)
        case Deleted      => BSONInteger(9)
      }
    }
  }

}
