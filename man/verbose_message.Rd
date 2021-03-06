% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index.R
\name{verbose_message}
\alias{verbose_message}
\title{Message function that can be globally muted}
\usage{
verbose_message(str, ...)
}
\description{
Prints a message if  \code{option("safedata.verbose")}  is set
to TRUE. Note that individual expressions can be muted using
\code{suppressMessages()} but this mutes them globally.
}
\keyword{internal}
