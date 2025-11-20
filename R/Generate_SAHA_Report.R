#' Generate an HTML Report from SAHA Analysis
#'
#' This function renders an HTML report using the SAHA analysis results.
#' It determines whether the provided output file path is absolute or relative,
#' resolves the final path accordingly, and uses `rmarkdown::render()` to
#' generate the report.
#'
#' @param SAHA An object containing SAHA analysis results.
#' @param auto Optional. A dataframe output from AutoAnnotate().
#' @param output_file A string specifying the output file name or full path for the report.
#'
#' @return The function creates an HTML report at the specified location. It does not return an R object.
#'
#' @details
#' - If `output_file` is an absolute path (e.g., "/home/user/report.html" or "C:/report.html"),
#'   it is used directly.
#' - If it is a relative path (e.g., "results/report.html"), it is resolved
#'   against the current working directory using `file.path(getwd(), output_file)`.
#' - The report template is an internal R Markdown file located within the `SAHAdata` package.
#'
#' @examples
#'
#' @export
Generate_SAHA_Report <- function(SAHA, auto = NULL, output_file) {
   
   # This regex checks for paths starting with '/' (Unix-like) or 'X:/' (Windows)
   is_absolute_path <- grepl("^(/|[A-Za-z]:/)", output_file)
   
   if (is_absolute_path) {
      # If it's already an absolute path, use it directly
      final_output_path <- normalizePath(output_file, mustWork = FALSE)
   } else {
      # If it's a relative path, combine it with the current working directory
      final_output_path <- normalizePath(file.path(getwd(), output_file), mustWork = FALSE)
   }
   
   # Render the .Rmd file
   rmarkdown::render(input = system.file("extdata", "report_html.Rmd", package = "SAHAdata"),
                     output_file = final_output_path)
}
