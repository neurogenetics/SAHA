Generate_SAHA_Report <- function(SAHA, auto=NULL, output_file) {
  # render the .Rmd file
  rmarkdown::render(input = system.file("extdata", "report_html.rmd", package = "SAHAdata"),
                    output_file = output_file)
}
