\name{marginalmodelplots-package}
\alias{marginalmodelplots-package}
\alias{marginalmodelplots}
\alias{marginalmodelplot}
\docType{package}
\title{Marginal Model Plots}
\description{Marginal mean model plots for linear (lm) and generalized linear models (glm).  
	Including tools for bandwidth exploration.  Built on the work from the alr3 and locfit package.}
\details{
\tabular{ll}{
Package: \tab marginalmodelplots\cr
Type: \tab Package\cr
Version: \tab 0.6\cr
Date: \tab 2007-05-15\cr
License: \tab GPL version 2 or newer.\cr
}
	The primary function of interest is \code{\link{mmplot}}, which creates a layout 
	of marginal mean model plots for a given model.
	The package also includes functions for exploring marginal model plots and 
	experimentation with bandwidths through the functions
	\code{\link{mmpslider}} and \code{\link{mmpdualslider}}.\\
	marginalmodelplots is based on the work of marginal model plots by Cook and Weisburg, 
	and extends the their work in from the package  \code{\link[alr3:alr3-package]{alr3}}. 
	Unlike \code{\link[alr3:alr3-package]{alr3}} the marginalmodelplots package only performs
	mean curve smooths.  For variance curve estimation use \code{\link[alr3:alr3-package]{alr3}}.
}
\author{
	Andrew Redd

	Maintainer: Andrew Redd <aredd at stat.tamu.edu>
}
\references{
	S. Weisberg (2005), \emph{Applied Linear Regression}, third edition, Wiley, Chapter 8
}
\keyword{iplot}
\keyword{smooth}
\keyword{regression}
\keyword{dynamic}
\seealso{
 \code{\link[alr3:alr3-package]{alr3}}
}

