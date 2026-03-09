\encoding{utf8}
\name{patterndependence}
\alias{patterndependence}
\alias{plot.pattern}
\alias{print.pattern}

\title{
Ordinal Pattern Dependence
}
\description{
Calculates the ordinal pattern coefficient and related values.
}

\usage{
patterndependence(tsx,tsy,d=3,block=FALSE,first=TRUE,
tiesmethod=c("random","first"),ordinalcor=c("standard","positive","negative"))

\method{plot}{pattern}(x, ...)
\method{print}{pattern}(x, ...)
}

\arguments{
	\item{tsx}{numeric vector representing the first univariate time series.}
	\item{tsy}{numeric vector representing the second univariate time series.}
	\item{d}{numeric value determining the length of the ordinal pattern.}
	\item{block}{logical value determining whether patterns are calculated on disjoint blocks or overlapping blocks.}
	\item{first}{logical value indicating which observartions are dropped if \code{block == TRUE} and the time series length is no multiple of \code{d}.}
	\item{tiesmethod}{character string specifying how ties, that is equal values, are treated, see `Details'.}
	\item{ordinalcor}{character string specifying which ordinal pattern coefficient is output, see `Details'.}
    \item{x}{object of class \code{"pattern"}, which is the output of \code{patterndependence}.}
    \item{...}{further arguments passed to the internal plotting function.}
}

\value{
Object of class \code{"pattern"} containing the following values:
	\item{patterncoef}{ordinal pattern coefficient.}
	\item{numbequal}{number of equal ordinal patterns.}
	\item{numbopposite}{number of opposite ordinal patterns.}
	\item{PatternXz}{number of ordinal patterns in first time series.}
	\item{PatternYz}{number or ordinal patterns in second time series.}
	\item{coding}{coding of the ordinal patterns, used in \code{PatternXz} and \code{PatternYz}.}
    \item{PatternX}{numeric vector representing the time series of patterns in \code{tsx}.}
    \item{PatternY}{numeric vector representing the time series of patterns in \code{tsy}.}
    \item{tsx}{numeric vector representing the first univariate time series.}
    \item{tsy}{numeric vector representing the second univariate time series.}
    \item{maxpat}{number representing the maximal pattern code.}
    \item{ordinalcor}{character string specifying the type of ordinal pattern coefficient.}
    \item{tiesmethod}{character string specifying how ties are treated.}
    \item{block}{logical value determining whether patterns are calculated on disjoint blocks or overlapping blocks.}
    \item{d}{length of the ordinal pattern.}
    \item{tablesame}{numeric vector representing the number of coinciding patterns, apportioned into different patterns.}
    \item{tableopposite}{numeric vector representing the number of reflected patterns, apportioned into different patterns.}
    \item{indexsame}{logic vector indicating whether patterns in both time series coincide.}
    \item{indexopposite}{logic vector indicating whether patterns in both time series are reflected.}
}

\details{
The standard ordinal pattern coefficient is a non-parametric and robust measure of dependence between two time series. It is based on ordinal patterns, which are defined as sequences of ranks of \code{d} subsequent observations. This sequences of subseqent observations can either move one observation per time or a whole block of \code{d} observations. The former is preferred since it uses more information. If one chooses the later, one has to decide whether the first or the last observations are removed in case that the time series length is no multiple of \code{d}. With regard to equal values within a window of consecutive observations (ties), the argument \code{tiesmethod} determines the approach for computing the respective ordinal patterns. The ``\code{first}'' method is in favor of increasing patterns, whereas the default ``\code{random}'' puts the equal values in random order.

Beside the default standard ordinal pattern coefficient, which range from -1 to 1, one can also look at the positive and negative ordinal pattern coefficient, which roughly measures whether there are unusual many identical or opposite patterns in the time series.  

The plot function draws both time series and shows the six most frequent coinciding pattern with counts on the right. At the bottom, the location of these coinciding patterns is visualized.
}

\references{

Schnurr, A. (2014): An ordinal pattern approach to detect and to model leverage effects and dependence structures between financial time series, \emph{Statistical Papers}, vol. 55, 919--931.

Schnurr, A., Dehling, H. (2017): Testing for Structural Breaks via Ordinal Pattern Dependence, \emph{Journal of the American Statistical Association}, vol. 112, 706--720.
}


\author{
Alexander \enc{Dürre}{Duerre}, 
Angelika Silbernagel
}


\examples{
set.seed(1066)
patternobj <- patterndependence(rnorm(100),rnorm(100))
plot(patternobj)
}
