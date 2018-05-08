package exports.results

import org.scalajs.dom.raw.HTMLSelectElement
import org.scalajs.jquery._

import scala.scalajs.js
import scala.scalajs.js.annotation.{ JSExport, JSExportTopLevel }

@JSExportTopLevel("TemplateAlignment")
class TemplateAlignment(tool: String, forwardingEnabled: Boolean = true) {

  private val templateAlignmentModal = jQuery("#templateAlignmentModal")
  private val textArea               = templateAlignmentModal.find("textarea.alignmentTemplateTextArea")
  private val toolSelect             = templateAlignmentModal.find(".alignmentTemplateToolSelect")

  if (forwardingEnabled) {
    templateAlignmentModal.find(".hide-for-forwarding-disabled").show()
    toolSelect.on(
      "click", { selectElement: HTMLSelectElement =>
        {
          val forwardData = textArea.value().toString
          val tool        = selectElement.value
          if (tool != "") {
            templateAlignmentModal.asInstanceOf[exports.facades.JQuery].foundation("close")
            Forwarding.simple(tool, forwardData)
          }
        }
      }: js.ThisFunction
    )
  } else {
    templateAlignmentModal.find(".hide-for-forwarding-disabled").hide()
  }

  templateAlignmentModal.on("open.zf.reveal", () => {
    toolSelect.value("")
    textArea.value("")
  })

  @JSExport
  def get(jobID: String, accession: String): Unit = {
    textArea.asInstanceOf[exports.facades.JQuery].LoadingOverlay("show")
    val acc = accession.replace("#", "%23")
    val extension = tool match {
      case "hhpred"  => "a3m"
      case "hhblits" => "ra3m"
      case "hhomp"   => "fas"
    }

    jQuery
      .ajax(
        js.Dictionary(
            "url" -> s"/results/templateAlignment/$jobID/$acc"
          )
          .asInstanceOf[JQueryAjaxSettings]
      )
      .done((_: js.Any, _: js.Any, _: JQueryXHR) => {
        jQuery
          .ajax(
            js.Dictionary(
                "url" -> s"/files/$jobID/$acc.$extension"
              )
              .asInstanceOf[JQueryAjaxSettings]
          )
          .done((data: js.Any, _: js.Any, _: JQueryXHR) => {
            textArea.value(data.toString)
          })
          .fail((jqXHR: JQueryXHR, textStatus: js.Any, errorThrow: js.Any) => {
            println(s"jqXHR=$jqXHR,text=$textStatus,err=$errorThrow")
            textArea.value("Sorry, failed to fetch Template Alignment.")
          })
          .always(() => {
            textArea.asInstanceOf[exports.facades.JQuery].LoadingOverlay("hide")
          })
      })
      .fail((jqXHR: JQueryXHR, textStatus: js.Any, errorThrow: js.Any) => {
        println(s"jqXHR=$jqXHR,text=$textStatus,err=$errorThrow")
        textArea.value("Sorry, failed to fetch Template Alignment.")
      })
      .always(() => {
        textArea.asInstanceOf[exports.facades.JQuery].LoadingOverlay("hide")
      })
  }
}
