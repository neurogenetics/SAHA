\name{Generate_SAHA_Report}
\alias{Generate_SAHA_Report}
\title{Generate SAHA HTML Report from Analysis Object}
\usage{
  Generate_SAHA_Report(SAHA, auto = NULL, output_file)
}
\description{
  This function renders an HTML report from a SAHA analysis object using a pre-defined R Markdown template. The report is saved to a user-specified path. If the path is relative, it is resolved against the current working directory.
}
\arguments{
  \item{SAHA}{An object containing SAHA analysis results.}
  \item{auto}{Optional. Reserved for future use to enable automatic mode configuration. Currently unused.}
  \item{output_file}{A string specifying the name or full path of the output HTML report. Relative paths will be resolved to the working directory.}
}
\details{
  This function checks whether the specified \code{output_file} is an absolute or relative path. It then uses \code{rmarkdown::render()} to generate a report using an internal template located in the \code{extdata} folder of the \code{SAHAdata} package.
}
\examples{
  # Save report to a file in the current directory
  Generate_SAHA_Report(SAHA = my_saha_object, output_file = "SAHA_Report.html")
  
  # Save report to an absolute path
  Generate_SAHA_Report(SAHA = my_saha_object, output_file = "/Users/username/Desktop/SAHA_Report.html")
}
