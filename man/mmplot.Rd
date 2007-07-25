\name{mmplot}
\alias{mmplot}
\alias{mmplot1}
\alias{mmplot1.lm}
\alias{mmplot1.glm}
\title{Marginal model plot}
\description{
Marginal model plots are a non parametric form of verifying model adequacy.
The Structure of a marginal model plot is a plot with the response variable
along the vertical axis, and a predictor variable, or any linear combination
of predictor variables allong the horizontal axis.  Two non parametric smooths
are added to the plots.  One for the response on the predictor and another for
the fitted values on the predictor.  If the two curves match closely then the 
model is determined to be adaquate.

The \code{mmplot} fuction, as implimented here, is capable of computing an optimal 
bandwidth using generalized cross validation (\code{gcv}), Akaike's information criteria (\code{aic}), 
or Akaike's information criteria corrected (\code{aicc}).
}
\usage{
mmplot(object,exclude=NULL,layout=NULL,ask,...)

mmplot1(object, ...)
mmplot1.lm(object,u=predict(object),
	label=deparse(substitute(u)), locfit.control=list(nn=2/3),
    colors=c("blue","red"), ...)
mmplot1.glm(object,u=predict(object),mean=TRUE,
    label=deparse(substitute(u)), locfit.control=list(nn=2/3),
    colors=c("blue","red"),family=NULL,link=NULL,yhat.autoquasi=TRUE,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A regression object, of type either lm or glm,
   for which there is a \code{predict}
  method defined. }
  \item{exclude}{A vector of indices of variables for which a plot is not
  required.  For example, \code{c(2,4)} would exclude the second and fourth
  terms in the model, excluding the intercept.}
  \item{layout}{A reasonable layout for the plots in the window is
  determined by the program.  If you don't like the default you can set your
  own layout:  \code{c(2,3)} means two rows and three columns.}
  \item{ask}{If TRUE, ask before clearing the graph window to draw more plots.}
  \item{u}{ The quantity to be plotted on the horizontal axis.  The
  default is the predicted values \code{predict(object)}.  Can be any other
  vector of length equal to the number of observations in the object. }
  \item{mean}{ If TRUE, compare mean smooths.  Depreciated, should allways be left as TRUE.}
  \item{label}{ label for horizontal axis }
  \item{locfit.control}{A list of control parameters to pass to the locfit function to adjust
			the nature of the curve.  the primary purpose of this will be to set the bandwidth 
			or nearest neighbor parameter for smoothing, although any of the \code{\link[locfit]{locfit}}
			parameters may be set here.  To use a optimally computed bandwidth use \code{h='aic'} or \code{nn='gcv'}
			for an AIC optimal fixed bandwidth or a GCV optimal nearest neighbor, respectivly.}
  \item{colors}{ colors for response and fitted values smooths, respectively }
  \item{family}{Overrides the family used for the local glm fit in glm models.  }
  \item{link}{Overrides the link used for the local glm fit in glm models. }
  \item{yhat.autoquasi}{Determines if the quasi version of the family should be used to fit the fitted values curve.  The default is True, since the fitted values do not have constant variance.}
  \item{\dots}{ additional arguments passed to \code{plot} }
}   
\details{
\code{mmplot1} draws one marginal model plot.  \code{mmplot} draws all marginal model plots,
versus each of the terms in the model and versus fitted values.  \code{mmplot}
should not be used if some of the terms in the model are factors or
interactions.

The \code{mmplot} implimentation here uses the locfit package to generate the non parametric smooths.
A \code{lm} model is equivalent to using a local polynomial smooth, but with \code{glm} models
the smooth is a local glm fit of the same family as the model, unless overridden by user.
}
\value{
\code{mmplot1} functions return (invisible) an object with the smothing parameters used.  
The first componenent being the parameters for the data (named \code{y}), the second for the fitted values (named \code{yhat}).  
The components of the smoothing parameters are first the bandwidth (\code{h}) if used, otherwise \code{NULL}. Second the nearest neighbor parameter(\code{nn}) if used, otherwise \code{NULL}.
The \code{mmplot} function returns an unnamed list of the smoothing parameters in the order the graph appear.
}
\note{
A large number of warnings can be produced from locfit when running mmplot, especially with optimized smoothing parameters.
Warnings that are passed from the locfit procedures are for the final fit.
}

\seealso{\code{\link[locfit]{locfit}}, \code{\link{plot}} }

\references{ S. Weisberg (2005), \emph{Applied
Linear Regression}, third edition, Wiley, Chapter 8}
\author{Andrew Redd, \email{aredd at stat.tamu.edu}}

\examples{
data(Bordeaux)
m<-lm(Price~ParkerPoints+CoatesPoints+FirstGrowth+Pomerol+VintageSuperstar+CultWine,data=Bordeaux)
mmplot(m,locfit.control=list(nn=2/3))

data(MichelinNY)
m4 <- glm(InMichelin~Food+Decor*Service+log(Cost),family=binomial(),data=MichelinNY)
mmplot(m4)
mmplot(m4,locfit.control=list(nn=.7))
mmplot(m4,locfit.control=list(nn='aic'))
} 
\keyword{hplot }% at least one, from doc/KEYWORDS
\keyword{regression}
