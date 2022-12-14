% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{generate_multilevel_data}
\alias{generate_multilevel_data}
\title{Functions used to generate data in the various scenarios in accompanying paper}
\usage{
generate_multilevel_data(
  scenario = 1,
  grid = seq(from = 0, to = 1, length.out = 48),
  N,
  p = rep(0.5, 2),
  J = 7,
  sigma = 0,
  binary = T,
  Ys = rbinom(N, 1, 0.5),
  return_scores_too = F
)
}
\arguments{
\item{scenario}{Which scenario to generate data for: 1=A, 2=B, 3=C In paper}

\item{grid}{which set of points within [0,1] will data be observed on}

\item{N}{Total number of subjects to generate data for}

\item{p}{vector showing distribution of the two groups}

\item{J}{number of realizations}

\item{sigma}{variance in the latent curves}

\item{binary}{Reports whether to return binary or latent curves}

\item{Ys}{Groups of the individuals}

\item{return_scores_too}{T/F values to determine whether to return score coefficients too}
}
\value{
: Binary or latent curves
}
\description{
Functions used to generate data in the various scenarios in accompanying paper
}
