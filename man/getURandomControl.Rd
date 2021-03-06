% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getUStart.R
\name{getURandomControl}
\alias{getURandomControl}
\title{Random starting weights}
\usage{
getURandomControl(nStart = 10, p = 0.8)
}
\arguments{
\item{nStart}{Number of random starts.}

\item{p}{The amount of weight to be assigned to a randomly selected class
for each observation. The remaining weight is split between the remainning classes.
\code{p} should be between zero and one, and is intended to be closer to one.}
}
\value{
An S3 object of class 'uRandomControl' containing \code{nStart} and \code{p}.
}
\description{
Constructor function for the list of parameters to be passed to the \code{uRandomControl}
argument of the \code{cac} function. For this starting method, a class is randomly
selected for each observation and given starting weight \code{p}. The remaining
weight is divided among the other classes.
}
