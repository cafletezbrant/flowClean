\name{clean}
\alias{clean}

\title{
  clean.  For cleaning flow cytometry data.
}
\description{
  This function uses compositional data analysis to identify errant
  collection events.
}
\usage{
clean(fF, vectMarkers, filePrefixWithDir, ext, binSize=0.01,
  nCellCutoff=500, announce=TRUE, cutoff="median", diagnostic=FALSE, fcMax=1.3)
}
\arguments{
  \item{fF}{ flowFrame object containing experimental data to be
    cleaned. }
  \item{vectMarkers}{ A vector of indices representing flow parameters
  to be examined.  These are considered as columns in the data matrix in
  which cells are rows and parameters are columns.  Generally this
  vector excludes indices for various `scatter` parameters
  (e.g. `FSC-A`)}
  \item{filePrefixWithDir}{ A string containing at least the desired
  name for the output flow file generated.  Can include directory
  structure and folder ('/' or '\') characters. }
  \item{ext}{ The file extension for the output flow file. }
  \item{binSize}{ A number in [0,1]; represents the fraction of duration
    of collection per bin.}
  \item{nCellCutoff}{ An integer; represents the minimum number of cells
    a population must have to be included in analysis.}
  \item{cutoff}{ Method for determining threshold for parameter.  Can be
  "median" (default) or in [0, 1], which is interpreted as a
  perecntile.  Integers > 1 will be interpreted as the fluorescence
  value to be used for a threshold.}
  \item{announce}{ Print completion messages. }
  \item{fcMax}{ Maximum allowable increase relative to presumed 'good' data. }
  \item{announce}{ If TRUE, will print message to screen if errors
    detected.}
  \item{diagnostic}{ If TRUE, will make PNG of populations in time
    bins, and save with same prefix as specified in filePrefixWithDir.}
  \item{returnVector}{ If desired, only return vector indicating if a
    given cell is 'good' or 'bad'.}
  \item{nstable}{ The number of stable populations required to be
  observed during the duration of an experiment.  Default is 5.}

}

\references{
  Fletez-Brant C, Spidlen J, Brinkman R, Roederer M and Chattopadhyay P.
  flowClean: Automated identification and removal of fluorescence
  anomalies in flow cytometry data.  Cytometry Part A, 2016.
}
\author{Kipper Fletez-Brant}

\seealso{The package vignette.}
\examples{
data(synPerturbed)
synPerturbed.c <- clean(synPerturbed, vectMarkers=c(5:17),
  filePrefixWithDir="sampleName", ext="fcs")
}
\keyword{quality control}
