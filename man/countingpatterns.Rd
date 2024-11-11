\encoding{utf8}
\name{countingpatterns}
\alias{countingpatterns}
\alias{print.patterncounts}

\title{
Empirical Ordinal Pattern Distribution
}
\description{
Calculates the empirical ordinal pattern distribution.
}

\usage{
countingpatterns(tsx,d=3,block=FALSE,first=TRUE,tiesmethod=c("random","first"),
generalized=FALSE)

\method{print}{patterncounts}(x, ...)
}

\arguments{
	\item{tsx}{numeric vector representing the univariate time series.}
	\item{d}{numeric value determining the length of the ordinal pattern.}
	\item{block}{logical value determining whether patterns are calculated on disjoint blocks or overlapping blocks.}
	\item{first}{logical value indicating which observartions are dropped if \code{block == TRUE} and the time series length is no multiple of \code{d}.}
	\item{tiesmethod}{character string specifying how ties, that is equal values, are treated if \code{generalized == FALSE}, see `Details'.}
	\item{generalized}{logical value determining whether classical ordinal patterns or their generalization with regard to ties are considered, see `Details'.}
  \item{x}{object of class \code{"patterncounts"}, which is the output of \code{countingpatterns}.}
  \item{...}{further arguments passed to the internal plotting function.}
}

\value{
Object of class \code{"patterncounts"} containing the following values:
	\item{patterncounts}{absolute frequencies of ordinal patterns.}
	\item{allpatterns}{list of all ordinal patterns considered.}
    \item{d}{length of the ordinal pattern.}
    \item{generalized}{logical value determining whether classical ordinal patterns or their generalization with regard to ties are considered.}
    \item{tiesmethod}{character string specifying how ties are treated.}
    \item{block}{logical value determining whether patterns are calculated on disjoint blocks or overlapping blocks.}    
}

\details{
Ordinal patterns, which are defined as sequences of ranks of \code{d} subsequent observations, are a useful tool to describe the dependence within or between time series. That sequences of subseqent observations can either move one observation per time or a whole block of \code{d} observations. The former is preferred since it uses more information. If one chooses the later, one has to decide whether the first or the last observations are removed in case that the time series length is no multiple of \code{d}. With regard to equal values within a window of consecutive observations (ties), the argument \code{tiesmethod} determines the approach for computing the respective ordinal patterns. The ``\code{first}'' method is in favor of increasing patterns, whereas the default ``\code{random}'' puts the equal values in random order.

Beside the classical ordinal patterns, one can also consider the generalized version proposed by Schnurr and Fischer (2022), where the ordinal information of stagnation in the case of ties is also included by taking into account a larger set of patterns.
}

\references{

Schnurr, A. (2014): An ordinal pattern approach to detect and to model leverage effects and dependence structures between financial time series, \emph{Statistical Papers}, vol. 55, 919--931.

Schnurr, A., Dehling, H. (2017): Testing for Structural Breaks via Ordinal Pattern Dependence, \emph{Journal of the American Statistical Association}, vol. 112, 706--720.

Schnurr, A., Fischer, S. (2022): Generalized ordinal patterns allowing for ties and their application in hydrology, \emph{Computational Statistics & Data Analysis}, vol 171, 107472.
}


\author{
Angelika Silbernagel
}


\examples{
set.seed(1066)
countingpatterns(rnorm(100))
countingpatterns(rpois(100,1), generalized=TRUE)
}
