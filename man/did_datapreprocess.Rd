% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperfuncs.R
\name{did_datapreprocess}
\alias{did_datapreprocess}
\title{Data cleaning and other helper functions
Julia Hatamyar, 11 May 2022}
\usage{
did_datapreprocess(outcome, group, time, id_name, data)
}
\arguments{
\item{outcome}{A vector of outcomes}

\item{group}{A vector indicating group membership}

\item{time}{A vector indicating time periods}

\item{id_name}{A vector indicating data ID (state, municipality, etc)}

\item{data}{A dataframe containing your data}
}
\value{
an intermediate processed step to be used within the MLDID function
}
\description{
Function to process arguments passed to the main methods
}
\details{
========================================
}
