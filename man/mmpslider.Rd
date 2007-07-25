\name{mmpslider}
\alias{mmpslider}
\alias{mmpslider.lm}
\alias{mmpslider.glm}
\alias{mmpdualslider}
\alias{mmpdualslider.lm}
\alias{mmpdualslider.glm}
\title{Margial Model Plots with Single and Dual Slider Bars}
\description{
\code{mmpslider} Creates an interactive marginal model plot for the supplied model and predictor variable.  
The graph provides an area at the bottom of the plot with a slider bar that can be used to interactivly adjust
the fixed bandwidth used in generating the nonparametric smooths in the plot.
\code{mmpdualslider} behaves behaves like \code{\link{mmpslider}} but creates two slider bars one for each fit, the response and the fitted valeues.  The slider bar is not actually a slider bar but an interactive region of the graph that will adjust the fixed bandwidth(s) of the non parametric smooth when clicked.  It does not actually slide.}
\usage{
mmpslider(m, ...)
mmpslider.lm(m,pred=predict(m), bw=NULL , label=deparse(substitute(u)), colors=c('blue','red'),...)
mmpslider.glm(m,pred=predict(m), bw=NULL , label=deparse(substitute(u)), colors=c('blue','red'), family=NULL , link=NULL, ...)
mmpdualslider(m, ...)
mmpdualslider.lm(m, pred = predict(m), bw = NULL, label = deparse(substitute(pred)), colors = c("blue", "red"),...)
mmpdualslider.glm(m, pred = predict(m), bw = NULL, label = deparse(substitute(pred)), colors = c("blue", "red"), family = NULL, link=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Model to use to generate marginal model plot from.}
  \item{pred}{The predictor variable or cobination of variables to use for the plot.  Defaults to the linear predictor of the model.}
  \item{bw}{Vector of length 3 that contains the minimun allowable bandwidth, the starting bandwidth and the maximum allowable bandwidth.  This setting sets up the range of the slider bar. The default for the minumum bandwidth is the maximum of the differences of sequential points.  The default starting value is the aic optimum value and te default maximum is twice the aic optimum bandwidth.  Either a vector of length 3 must be provided of left null.  \code{mmpdualslider} uses a local glm fit for marginal model plots with glm models.}
  \item{label}{The x axis label.  So named to reflect the same parameter in similar mmp code.}
  \item{colors}{Vector of length 2 containg the colors for the smooth against the response and the fitted values, respectivly.}
  \item{family}{Overrides the family used for the local glm fit in glm models.  }
  \item{link}{Overrides the link used for the local glm fit in glm models. }
  \item{\dots}{Passes arguments onto other functions}
}
\details{
To use \code{mmpdualslider} run the function then click in the slider bar area forto adjust the bandwidth
 of that respective curve. The slider bars are color coordinated with the curve that they represent as well
 as being labeled in the slider indicator.  Clicking the stop button, or outside of the slider bars area will
 finalize the plot.\\
The plot is a custom plot generated from scratch with the grid utilities.  Any further customization of the 
plot will have to be used with in the grid framework.
}
\value{none. Used for side effects of generating a plot.}
\references{
	S. Weisberg (2005), \emph{Applied
	Linear Regression}, third edition, Wiley, Chapter 8
}
\author{Andrew Redd <aredd at tamu.edu>}

\seealso{\code{\link{marginalmodelplots-package}}}
\examples{
\dontrun{
data(Bordeaux)
names(Bordeaux)
attach(Bordeaux)
m<-lm(Price~ParkerPoints+CoatesPoints+FirstGrowth+Pomerol+VintageSuperstar+CultWine,data=Bordeaux)
mmpdualslider(m,ParkerPoints)
detach(Bordeaux)

data(MichelinNY)
attach(MichelinNY)
m <- glm(InMichelin~Food+Decor+Service+Cost+log(Cost)+Food:Decor,family=binomial(),data=MichelinNY)
mmpdualslider(m)
mmpdualslider(m,Food,bw=c(1,3,10))
mmpdualslider(m,Decor,bw=c(1,3,10))
detach(MichelinNY)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{iplot}
\keyword{smooth}
\keyword{regression}
\keyword{dynamic}
