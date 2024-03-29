% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate-lambda.R
\name{coxintensity2matrix}
\alias{coxintensity2matrix}
\title{Evaluate lambda}
\usage{
coxintensity2matrix(x, nx = 2^8, ny = NULL, bbox, at = NULL, W, ...)
}
\arguments{
\item{x}{object from 'lambda'-function}

\item{nx}{if 'at' missing, make a grid with this many x-steps (only in 2d!)}

\item{ny}{like nx but for y-dimension (can be omitted)}

\item{bbox}{bounding box for the window in case grid used}

\item{at}{where to evaluate, matrix of coordinates}

\item{W}{alternative to bbox, owin-object (2D!)}

\item{...}{omitted.}
}
\description{
Evaluate lambda
}
