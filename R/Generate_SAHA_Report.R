#' Generate a Summary HTML Report for a SAHA Object
#'
#' This function renders a pre-built RMarkdown report using the `report_html.Rmd` template bundled within the SAHAdata package. The report summarizes key analysis outputs stored  in the provided SAHA object.
#'
#' @param SAHA A SAHA object containing analysis results (not directly passed to the RMarkdown file,but typically expected to be loaded in the report environment).
#' @param auto Optional. A placeholder for potential future use, e.g., to toggle auto-annotation output.
#' 
#' @param output_file A file path (including filename) where the rendered HTML report should be saved.
#'
#' @return No return value. The function generates and saves an HTML file to `output_file`.
#'
#' @details
#' The report template is located in the "extdata" directory of the `SAHAdata` package. 
#' Users can customize this template if needed by copying and modifying it separately.
#'
#' @examples
#' Generate_SAHA_Report(SAHA = my_saha_object, output_file = "SAHA_summary.html")
#'
#' @export
Generate_SAHA_Report <- function(SAHA, auto=NULL, output_file) {
  # render the .Rmd file
  rmarkdown::render(input = system.file("extdata", "report_html.Rmd", package = "SAHAdata"),
                    output_file = output_file)
}
