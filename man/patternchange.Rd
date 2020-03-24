\encoding{utf8}
\name{patternchange}
\alias{patternchange}
\alias{plot.change}

\title{
Changepoint Detection Using Ordinal Patterns
}
\description{
Test for a change in the dependence structure of two time series using ordinal patterns
}

\usage{
patternchange(tsx,tsy,h=2,conf.level,weight=TRUE,weightfun=NULL,bn=log(length(tsx)),
    kernel=function(x){return(max(0,1-abs(x)))})

\method{plot}{change}(x, ...)
}

\arguments{
	\item{tsx}{numeric vector of first univariate time series.}
	\item{tsy}{numeric vector of second univariate time series.}
	\item{h}{numeric value determining the length of ordinal pattern.}
    \item{conf.level}{numerical value indicating the confidence level of the test.}
	\item{weight}{logical value indicating whether one uses weights of the L1 norm or the empirical probability of identical patterns; see details.}
	\item{weightfun}{function which defines the weights given the L1 norm between the patterns if \code{weight=TRUE}. If no weight-function is given, the canonical weight function is used; see details.}
	\item{bn}{numerical value determining the bandwidth of the kernel estimator used to estimate the long run variance.}
	\item{kernel}{kernel function for estimating the long run variance.}
	\item{x}{object of class \code{"change"}}
	\item{...}{further arguments passed to the internal plotting function (\code{plot}).}
}


\value{
Object with classes \code{"change"} and \code{"htest"} containing the following values:
	\item{statistic}{the value of the test statistic. Under the null the test statistic follows asymptotically a Kolmogorov Smirnov distribution.}
	\item{p.value}{the p-value of the test.}
	\item{estimate}{the estimated time of change.}
	\item{null.value}{the jump height of the at most one change point model, which is under the null hypothesis always 0.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a characters string describing the test.}
	\item{trajectory}{the cumulative sum on which the tests are based on. Could be used for additional plots.}
}

\details{
Given two timeseries \code{tsx} and \code{tsy} a cusum type statistic tests whether there is a change in the patter dependence or not. The test is based on a comparison of patterns of length \code{h+1} in \code{tsx} and \code{tsy}. One can either choose the number of identical patterns (\code{weight=FALSE}) or a metric that is defined by the \code{weightfun} argument to measure the difference between patterns (\code{weight=TRUE}). If no (\code{weightfun}) is given, the canonical weightfunction is used, which equals 1 if patterns are identical and 0 if the L1 norm of their difference attains the maximal possible value. The value is linear interpolated in between.

The procedure depends on an estimate of the long run variance. Here a kernel estimator is used. A kernel function and a bandwidth can be set using the arguments \code{kernel} and \code{bn}. If none of them is given, the bartlett kernel with a bandwidth of \code{log(n)}, where \code{n} equals the length of the timeseries, is used.  
}

\references{

Schnurr, A. (2014): An ordinal pattern approach to detect and to model leverage effects and dependence structures between financial time series, \emph{Statistical Papers}, vol. 55, 919--931.

Schnurr, A., Dehling, H. (2017): Testing for Structural Breaks via Ordinal Pattern Dependence, \emph{Journal of the American Statistical Association}, vol. 112, 706--720.
}

\seealso{
Estimation of the pattern dependence is provided by \code{\link{patterndependence}}.
}

\author{
Alexander \enc{Dürre}{Duerre}
}


\examples{
set.seed(1066)
a1 <- cbind(rnorm(100),rnorm(100))
a2 <- rmvnorm(100,sigma=matrix(c(1,0.8,0.8,1),ncol=2))
A <- rbind(a1,a2)
testresult <- patternchange(A[,1],A[,2])
plot(testresult)
testresult
}
