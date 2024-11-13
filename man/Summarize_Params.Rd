\name{Summarize_Params}
\alias{Summarize_Params}
\title{Collect Parameter Settings for SAHA Analysis}
\usage{
   Summarize_Params(ann)
}

\description{
   The `Summarize_Params` function gathers and organizes parameter settings used in various SAHA analyses, such as marker-based filtering and self-similarity, and stores them in a structured data frame within the `ann@params$summary` slot. This centralizes all analysis parameters for easier tracking and reproducibility.
}

\arguments{
   \item{ann}{A SAHA object that contains analysis parameter settings within the `ann@params` slot. The function updates `ann@params$summary` with a structured data frame of these settings.}
}

\value{
   Returns the input `ann` object with a new or updated data frame in `ann@results$params`, listing parameter settings for SAHA analyses.
}

\examples{
   # Collect SAHA parameters and store in ann@params$summary
   ann <- Summarize_Params(ann)
}
