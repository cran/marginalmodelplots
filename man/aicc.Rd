\name{aicc}
\alias{aicc}
\title{Akaike's Information Criteria - Corrected for Non-Parametic Models}
\description{
Calculates Akaike's information criteria corrected for non parametric models made with \code{\link[locfit]{locfit}}.
}
\usage{
aicc(x,\dots,pen=2)
}
\arguments{
  \item{x}{ A model or model formula.}
  \item{\dots}{ Other arguments to locfit}
  \item{pen}{ Penalty for the degrees of freedom term.}
}   
\value{
vector of length 4 with the value of the aicc in the fourth spot.
}

\seealso{\link[locfit]{aic}, \link[locfit]{gcv}}

\references{
Chiou, J. and Tsai, C. (2006). Smoothing parameter selection in quasi-likelihood models.
\emph{Nonparametric Statistics}, 18(3), 307-314.
}
\author{Andrew Redd,  \email{aredd at stat.tamu.edu}}

\keyword{models}
\keyword{regression}
