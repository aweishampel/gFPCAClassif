% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapping_functions.R
\name{gmFPCA_predict}
\alias{gmFPCA_predict}
\title{Function for predicting the groups for new gmFPCA Classifier}
\usage{
gmFPCA_predict(
  gmFPCA.model,
  X_dat_m_new,
  covariates_new = NA,
  gar_covariates_new = NA
)
}
\arguments{
\item{gmFPCA.model}{Trained gmFPCA}

\item{X_dat_m_new}{N_new x m matrix of binary data}

\item{covariates_new}{N_new x Q data frame of covariates}

\item{gar_covariates_new}{N_new X Q_g data frame of covariates for gAR model}
}
\value{
predicted group values for the N_new users
}
\description{
Function for predicting the groups for new gmFPCA Classifier
}
