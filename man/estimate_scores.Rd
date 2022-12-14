% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{estimate_scores}
\alias{estimate_scores}
\title{Solve the GLMM based on the Estimates of the scores for the multilevel scenario}
\usage{
estimate_scores(
  data_formatted,
  s_mat,
  J,
  I,
  eigen_vals1,
  eigen_vals2,
  eigen_funcs1,
  eigen_funcs2,
  mu_t
)
}
\arguments{
\item{data_formatted}{NJ x D binary matrix}

\item{s_mat}{Matrix of estimated latent state values S}

\item{J}{Number of curves per subject}

\item{I}{N Number of subjects}

\item{eigen_vals1}{Vector  of Eigenvalues of the first level (individual) effects}

\item{eigen_vals2}{Vector of Eigenvalues of the second level (daily) effects}

\item{eigen_funcs1}{Matrix K1 x D of eigenfunctions of first level}

\item{eigen_funcs2}{Matrix K2 x D of eigenfunctions of second level}

\item{mu_t}{Mean function vector}
}
\value{
First Level Scores for every individual
}
\description{
Solve the GLMM based on the Estimates of the scores for the multilevel scenario
}
