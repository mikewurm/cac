% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getUStart.R
\name{getUSubsetControl}
\alias{getUSubsetControl}
\title{Random subset starting values}
\usage{
getUSubsetControl(nStart = 10, npp = NULL, maxiter = 100)
}
\arguments{
\item{nStart}{Number of random starts.}

\item{npp}{Number of observations randomly selecting for each starting subset.}

\item{maxiter}{Maximum attempts before resorting to getUStartRandom. It is possible
for the "subset" starting method to fail, for example if one of the subsets
contains only censored observations. If "subset" starts faill too frequently,
it may help to increase \code{npp}.}
}
\value{
An S3 object of class 'uSubsetControl' containing \code{nStart},
\code{npp}, and \code{maxiter}.
}
\description{
Constructor function for the list fo parameters to be passed to the \code{uSubsetControl}
argument of the \code{cac} function. This method randomly selects a subset of
observations of size \code{npp} for each class. Beta for each class is initialized
to the Cox model fit within the corresponding subset. Observations are drawn
without replacement. The mixture probabilities are initialized to be equal, and
the baseline hazard functions are each initialized to be the
Nelson-Aalen estimate from the full-data.
}
