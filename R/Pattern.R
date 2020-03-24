patternseq <- function(timeseries,h=2,block=FALSE,first=TRUE) {
h <- h+1
n <- length(timeseries)                         
if (block==FALSE) {                 # using point by point method
	numb <- n-h+1                        # how many patterns
	result <- matrix(0,ncol=h,nrow=numb)  
	pattern <- .Call("vergleich",timeseries,result,rep(0,h))
	return(h-pattern-1)
	}
if (block==TRUE) {                 # using disjoint blocks
	numb <- n%/%h                            
	remov <- n%%h                  
	if (first==TRUE)   {
		timeseries <- timeseries[(remov+1):n]}    # remove first values
	if (first==FALSE)  {
		timeseries <- timeseries[1:(n-remov)]}    # remove last values
	pattern <- matrix(nrow=numb,ncol=h)  
	for (i in 1:numb)   {
		pattern[i,] <- h-rank(timeseries[((i-1)*h+1):(i*h)],ties.method="random") 
	    	}
    	return(pattern)
	}
}


patterndependence <- function(tsx,tsy,h=2,block=FALSE,first=TRUE) {
h <- h+1

# compute pattern
PatternX <- patternseq(tsx,h=h-1,block=block,first=first)+1  # ranks
PatternY <- patternseq(tsy,h=h-1,block=block,first=first)+1  # ranks
PatternX <- as.numeric(PatternX%*%c(0,h^(0:(h-2))))   # transformation to one-dimensional object
PatternY <- as.numeric(PatternY%*%c(0,h^(0:(h-2))))

# count pattern 
indexsame <- PatternX==PatternY
numbsame <- sum(indexsame)  
    
maxPattern <- as.numeric(t(1:h)%*%c(0,h^(0:(h-2))))   
minPattern <- as.numeric(t(h:1)%*%c(0,h^(0:(h-2))))
PatternX <- PatternX-minPattern
PatternY <- PatternY-minPattern
maxPattern <- maxPattern-minPattern
indexopposite <- PatternX+PatternY==maxPattern
numbopposite <- sum(indexopposite)  # value of opposite patterns equal maxPattern! 
 
allpattern <- permutations(n=h,r=h)
codepattern <- as.numeric(allpattern%*%c(0,h^(0:(h-2))))-minPattern # value of all possible patterns
PatternXz <- factor(PatternX,levels=codepattern)    # transformation to factors
PatternYz <- factor(PatternY,levels=codepattern)

tablesame <- table(PatternXz[indexsame])
tableopposite <- table(PatternXz[indexopposite])

PatternXz <- table(PatternXz)                      # counting the patterns
names(PatternXz) <- 0:(length(codepattern)-1)
PatternYz <- table(PatternYz)
names(PatternYz) <- 0:(length(codepattern)-1)
Patternnamen <- 0:(length(codepattern)-1)
coding <- cbind(allpattern-1,codepattern)
coding <- coding[order(coding[,4]),]



# calculate coefficients
n <- sum(PatternXz)
q <- sum(PatternXz*PatternYz)/n^2
p <- numbsame/n
r <- numbopposite/n
s <- sum(PatternXz*rev(PatternYz))/n^2
alpha <- p-q
beta <- r-s
ordercor <- max(alpha/(1-q),0)-max(beta/(1-s),0)
result <- list(ordercor,alpha,beta,numbsame,numbopposite,PatternXz,PatternYz,coding,PatternX,PatternY,tsx,tsy,maxPattern,block,h,tablesame,tableopposite,indexsame,indexopposite)
names(result) <- c("patterncoef","alpha","beta","numbequal","numboppos","patterncounttsx","patterncounttsy","coding",
"patternseqtsx","patternseqtsy","tsx","tsy","maxpat","block","h","tablesame","tableopposite","indexsame","indexopposite")
class(result) <- "pattern"
return(result)
}

patternchange <- function(tsx,tsy,h=2,conf.level=0.95,weight=TRUE,weightfun=NULL,bn=log(length(tsx)),kernel=function(x) {return(max(0,1-abs(x)))}) {
h <- h+1

# Defining standard weight function
if(weight==TRUE) {
	if(is.null(weightfun)) {
		maxdif <- floor(h/2)*(floor(h/2)+1)+floor((h-1)/2)*(floor((h-1)/2)+1)
		weightfun <- function(x) return((maxdif-x)/maxdif)
		}
	}

# Computing Patterns
PatternX <- patternseq(tsx,h=h-1,block=FALSE)
PatternY <- patternseq(tsy,h=h-1,block=FALSE)

# Calculating transformed observations
L1norm <- apply(abs(PatternX-PatternY),1,sum)

if (weight==TRUE) {
	obs <- weightfun(L1norm)
	}	else{
	obs <- L1norm==0
	}

# Calculation of long-run-variance 
n <- length(obs)
weightv <- kernel(0:floor(bn)/bn)
acfv <- acf(obs,lag.max=floor(bn),type="covariance",plot=FALSE)$acf
sigma <- acfv[1]+2*sum(acfv[2:floor(bn)])

# Calculation of Cusum statistic

Tn <- 1/sqrt(n)*abs(cumsum(obs-mean(obs)))/sqrt(sigma)
changepoint <- which.max(Tn)
Tnmax <- max(Tn)
pvalue <- 1-pKS2(Tnmax)
cval <- c(-1,+1)*qKS(conf.level)
names(changepoint) <- "estimated change point"
names(Tnmax) <- "test statistic"
null.value <- 0
names(null.value) <- "possible level shift"
result <- list(p.value=pvalue,
statistic=Tnmax,
null.value=0,
critical.value = cval,
alternsative="two.sided",
method = paste("Test for a change in the pattern dependence"),
estimate=changepoint,
trajectory=Tn)
class(result) <- c("htest","change")
attr(result$critical.value, "conf.level") <- conf.level
return(result)
}

pKS2 <- function(x) {
if (x<=0.05) return(0)
p <- 1-2*sum((-1)^(2:1001)*exp(-2*(1:1000)^2*x^2))
return(p)
}


qKS <- function(p) {
  pKSm <- function(y) return(pKS2(y)-p)
  quan <- uniroot(pKSm, lower = 0.2, upper = 3)
  res <- quan$root
  return(res)
}

plot.change <- function(x, ...){
  plot(x$trajectory, main = x$method, type = "l", xlab = "Time", ylab = "Test statistic", ylim = range(x$trajectory, x$critical.value, na.rm = TRUE), ..., sub = paste("(estimated change point: ", x$estimate, "; pvalue: ", round(x$p.value,3), ")", sep=""))
  abline(v = x$estimate, lty = "solid", col = "red")
 	abline(h = x$critical.value, lty = "dashed", col = "blue")
}

plot.pattern <- function(x, ...) {
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))

layout(matrix(c(as.numeric(rbind(1,2:13)),14:25),ncol=2,byrow=TRUE),widths=c(0.8,0.1),heights=rep(1,18))
par(mar=c(1,1,0.5,0),xaxt="n",yaxt="n")
tsx <- x$tsx
tsy <- x$tsy
PatternX <- x$patternseqtsx
PatternY <- x$patternseqtsy
h <- x$h
n <- length(tsx)
rangv <- c(min(tsx,tsy),max(tsx,tsy))
plot(tsx,type="l",ylim=rangv, ...)
lines(tsy,col="blue")


anz <- x$tablesame
anz <- sort(anz,decreasing=TRUE)

for (j in 1:6) {

    indexplo <- (PatternX==as.numeric(names(anz[j])))&x$indexsame
    indexcoding <- as.numeric(names(anz[j]))
    indexcoding <- which(x$coding[,(h+1)]==indexcoding)
    
    par(mar=c(0.5,1,0.5,1),bty="o")
    plot(x$coding[indexcoding,1:h],type="l",xlab="",ylab="",main="",lwd=3)

    par(mar=c(0.5,1,0.5,1),bty="n")
    plot(1,0,type="n",main="",xlab="",ylab="")
    text(1,0,anz[j],col="darkgreen")
}


for (j in 1:6) {

    par(mar=c(0.5,1,0.5,0),bty="o")
    plot(1:n,rep(0,n),type="n",main="",xlab="",ylab="",ylim=c(0,1))
    
    indexplo <- (PatternX==as.numeric(names(anz[j])))&x$indexsame

    index <- which(indexplo==1)
    nn <- length(index)
    if (nn > 0) {
        if (x$block==FALSE) {
            for (i in 1:nn) lines(c(index[i],index[i]),c(0,1),col="darkgreen",lwd=2)
            }
        if (x$block==TRUE) {
            if (n < 1000) {
                for (i in 1:nn) polygon(c((index[i]-1)*h+1,index[i]*h,index[i]*h,(index[i]-1)*h+1),c(0,0,1,1),col="green",border=NA)
                }
            if (n>=1000) {
                for (i in 1:nn) lines(c(index[i]*h,index[i]*h),c(0,1),col="darkgreen",lwd=2)
                }
            }
        }

    indexcoding <- as.numeric(names(anz[j]))
    indexcoding <- which(x$coding[,h+1]==indexcoding)

    par(mar=c(0.5,1,0.5,1),bty="o")
    plot(x$coding[indexcoding,1:h],type="l",xlab="",ylab="",main="",lwd=3)
}
}

print.pattern <- function(x, ...) {
cat("\n")
cat(" standardized ordinal pattern coefficient: ",x$patterncoef,sep="")
cat("\n")
cat(" using ",x$h+1," consecutive observations and ",sep="")
if (x$block==TRUE) cat("nonoverlapping blocks")
if (x$block==FALSE) cat("overlapping blocks")
cat("\n")
cat("\n")
cat(" there are",x$numbequal,"coinciding and",x$numboppos,"reflected patterns")
cat("\n")
cat("\n")
spacen1 <- 2*x$h+3-8
cat(" pattern")
for(i in 1:spacen1) cat(" ")
cat("| coinciding | reflected")
cat("\n")
patterncoin <- x$tablesame
indexcoin <- order(as.numeric(names(patterncoin)))
patterncoin <- patterncoin[indexcoin]
patternrefl <- x$tableopposite
indexrefl <- order(as.numeric(names(patternrefl)))
patternrefl <- patternrefl[indexrefl]
spacen2 <- max(patterncoin)
spacen2 <- 10-nchar(spacen2)
spacen3 <- max(patternrefl)
spacen3 <- 9-nchar(spacen3)
n <- length(x$tablesame)
for (i in 1:n) {
cat(" (")
for (j in 1:x$h) {
cat(x$coding[i,j])
if (j<x$h) cat(",") else {
cat(") | ")
}
}
for (j in 1:spacen2) cat(" ")
cat(patterncoin[i])
cat(" | ")
for (j in 1:spacen3) cat(" ")
cat(patternrefl[i])
cat("\n")
}
cat("\n")
}




