% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hs.R
\name{CCfoldHS.calc}
\alias{CCfoldHS.calc}
\title{CCfoldHS.calc}
\usage{
CCfoldHS.calc(
  p.dir = getwd(),
  s.dir,
  c.dir,
  coln = "NormHpMChr",
  readfilename = "Hotspot.Table",
  constant = 0
)
}
\arguments{
\item{p.dir}{Parent directory}

\item{s.dir}{sample directory containing hotspot tables to be x in log2(x/y)}

\item{c.dir}{sample directory containing hotspot tables to be y in log2(x/y)}

\item{coln}{name of the column to be used in log2(x/y)}

\item{readfilename}{name of the read pattern}
}
\description{
Performs log2(x/y) for a column in a hotspot table
}
\examples{
CCfoldHS.calc(c.dir="/Users/georgebrown/Desktop/foldhsdev/oligosae2d",
s.dir="/Users/georgebrown/Desktop/foldhsdev/s1",
p.dir=getwd(),coln="NormHpM")
}
