% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{nb_updated_scores_only}
\alias{nb_updated_scores_only}
\title{Function: Grid search to estimate predicted values and estimate the values of h}
\usage{
nb_updated_scores_only(scores, classes, prior_g, scores_test, h = 1.06)
}
\arguments{
\item{scores}{N x K matrix of scores in the training set}

\item{classes}{Group labels vector of length N}

\item{prior_g}{vector of prior probability of being in each group, sums up to 1}

\item{scores_test}{N_test x K matrix of scores in the testing set}

\item{h}{multiplier for kernel based density function}
}
\value{
predicted classes for the Bayes Classifier
}
\description{
Function: Grid search to estimate predicted values and estimate the values of h
}
