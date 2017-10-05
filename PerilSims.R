### Perils of NHT
### R scripts
### October 18, 2016
###

## Loading and installing (if not already installed)
##  packages

packages<-function(x, repos="http://cran.r-project.org", ...){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x, repos=repos, ...)
    require(x,character.only=TRUE)
  }
}

base <- getwd()
dataDIR <- paste(base, "Data", sep="/")
plotDIR <- paste(base, "Figures", sep="/")

packages(arm)
packages(lattice)
packages(tikzDevice)
packages(rv)
packages(car)
## Simulation of the ``nonparametric'' change point model of Qian et al (2003)
## With a normal response and equal variance


chngp <- function(infile)
{ ## infile is a data frame with two columns
  ##  Y and X
    temp <- na.omit(infile)
    yy <- temp$Y
    xx <- temp$X
    mx <- sort(unique(xx))
    m <- length(mx)
    vi <- numeric()
    vi [m] <- sum((yy - mean(yy))^2)
    for(i in 1:(m-1))
            vi[i] <- sum((yy[xx <= mx[i]] - mean(yy[xx <=
                mx[i]]))^2) + sum((yy[xx > mx[i]] - mean(
                yy[xx > mx[i]]))^2)
    thr <- mean(mx[vi == min(vi)])
    return(thr)
}

## Using apply:

cart.sim <- function(n=30, n.sims=10000){
    simdata <- cbind(matrix(runif(n*n.sims), ncol=n),
                     matrix(rnorm(n*n.sims), ncol=n))
    return(mean(apply(simdata, 1, function(x){
        n <- length(x)
        temp <- data.frame(X=x[1:(n/2)], Y=x[(n/2+1):n])
        split <- chngp(temp)
        if (split==min(temp$X) | split==max(temp$X)) p.value=0.5
        else p.value <- t.test(Y~I(X<split), data=temp, var.equal=T)$p.value
        return(p.value<0.05)
    })))
}

typeIErr <- 0
j <- 0
n <- seq(20, 200, 20)
set.seed(123)
for (i in n){
    j <- j+1
    print(i)
    typeIErr[j] <- cart.sim(n=i)
}


inverts <- read.csv(paste(dataDIR, "InvertebratesGeomeansUTP6mo.csv", sep="/"))
xyplot(BCD~GM.UTP|DATE, data=inverts)

cart.sim2 <- function(x=inverts$GM.UTP, y=inverts$BCD, n.sims=10000){
    n <- length(x)
    reject <- 0
    for (i in 1:n.sims){
        split <- chngp(Df <- data.frame(X=sample(x,n,replace=F), Y=y))
        if (split==min(Df$X) | split==max(Df$X)) p.value=0.5
        else p.value <- t.test(Y~I(X<split), data=Df, var.equal=T)$p.value
        reject <- reject + (p.value<0.05)/n.sims
    }
    return(reject)
}

set.seed(123)
typeIBCD <- 0
j <- 0
for (i in levels(inverts$DATE)[-1]){
    temp <- inverts$DATE==i
    j <- j+1
    print(c(i, sum(temp),
            round(typeIBCD[j] <- cart.sim2(x=inverts$GM.UTP[temp],
                                                         y=inverts$BCD[temp]), 3)))
}

tol <- read.csv(paste(dataDIR, "macroinv.tolerant.csv", sep="/"))
typeITol <- 0
j <- 0
set.seed(123)
for (i in levels(tol$Date)){
    temp <- tol$Date==i
    j <- j+1
    print(c(i, sum(temp),
            round(typeITol[j] <- cart.sim2(x=tol$TP_6mo_gm[temp],
                                           y=logit(tol$tolerant[temp]/tol$Total[temp])), 3)))
}

sz <- c(table(inverts$DATE), table(tol$Date))[-1]
typeIErr2 <- c(typeIBCD, typeITol)

tikz(file=paste(plotDIR, "typeIfig1.tex", sep="/"),
     height=3, width=4.5, standAlone=F)
par(mar=c(3,3,0.5,0.125), mgp=c(1.75,0.125,0), las=1, tck=0.01)
plot(n, typeIErr, log="x", xlab="Sample Size", ylab="Type I Error Probability")
points(jitter(as.numeric(sz)), typeIErr2, pch=3)
dev.off()


### Function for calculating IVs along a gradient
###  Corrected a coding error in B&K for handling ties in
###  gradient values

T1 <- function(grad, spAbun, minsplit=5, Log=F){
## calculating IVs using species abundance along a gradient
    if (Log) spAbun <- log10(spAbun+1)
    ## in some cases, B&K used log10(abundance + 1)
    oo <- order(grad)
    spAbun <- spAbun[oo] ## sorting
    grad <- grad[oo]
    splits <- unique(grad)
    n <- length(splits)
    nsp <- n-2*minsplit+1 ## number of potential splits
    b <- as.numeric(spAbun>0)
    IV <- numeric()
    for (i in 1:nsp){
        x1 <- sum(spAbun[grad<=splits[minsplit+i-1]])
        x2 <- sum(spAbun[grad>=splits[minsplit+i]])
        x1 <- x1/(minsplit+i-1)
        x2 <- x2/(n-minsplit-i+1)
        tabun <- x1+x2
        if(tabun==0)tabun <- 1
        Ra1 <- x1/(tabun)
        Ra2 <- 1-Ra1
        Rf1 <- sum(b[grad<=splits[minsplit+i-1]])
        Rf2 <- sum(b[grad>=splits[minsplit+i]])
        Rf1 <- Rf1/length(b[grad<=splits[minsplit+i-1]])
        Rf2 <- Rf2/length(b[grad>=splits[minsplit+i]])
        IV[i] <- 100*max(c(Ra1*Rf1, Ra2*Rf2))
    }
    new.grad <- splits[minsplit:(n-minsplit)]
    split <- seq(1,length(new.grad))[IV==max(IV)]
    if (min(split)==1 | max(split)==length(new.grad)) no.chngp <- 1
    else no.chngp <- 0
    if (max(split)==length(new.grad))
        split<- new.grad[max(split)]
    else
        split <- mean(new.grad[c(split, split+1)])
    return (list(IVs=data.frame(splits=splits[minsplit:(n-minsplit)],
                     IndVal=IV),
                 split=c(split, no.chngp)))
}

### Function for permutation test for a specific split
### See notes above
perm <- function(numperm = 250, grad, spAbun, minsplit=5, Log=FALSE){
  oo <- order(grad)
  grad<- grad[oo]
  spAbun <- spAbun[oo]
  b <- as.numeric(spAbun >0)
  if (Log) spAbun <- log10(spAbun+1)
  split <- T1(grad, spAbun, minsplit, Log=Log)
  split.obs<-max(split$IVs$IndVal)
  cp <- split$split[1]
  n1 <- sum(grad<=cp)
  n2 <- sum(grad> cp)
  split.perm <- numeric()
  grad.count <- table(grad)
  grad.order <- as.numeric(ordered(grad))
  for (i in 1:numperm){
    perms <- sample(rep(1:2, c(n1, n2)), size=max(grad.order))
    while (sum(perms==1)==0 | sum(perms==2)==0)
      perms <- sample(rep(1:2, c(n1, n2)), size=max(grad.order))
    grad.od <- rep(perms, grad.count)
    i1=which(grad.od==1)
    i2=which(grad.od==2)
    x1 <- sum(spAbun[i1])
    x2 <- sum(spAbun[i2])
    x1 <- x1/length(unique(grad[i1]))
    x2 <- x2/length(unique(grad[i2]))
    tatl <- x1+x2
    if (tatl == 0) tatl<-1
    Ra1 <- x1/(tatl)
    Ra2 <- 1-Ra1
    Rf1 <- sum(b[i1])
    Rf2 <- sum(b[i2])
    Rf1 <- Rf1/length(i1)
    Rf2 <- Rf2/length(i2)
    split.perm[i] <- 100*max(c(Ra1*Rf1, Ra2*Rf2))
  }
  return(c(p.value=mean(split.perm >= split.obs),
              z=(split.obs-mean(split.perm))/sd(split.perm),
              perm.mu = mean(split.perm), perm.sd = sd(split.perm)))
}

## Type I error probability
typeIsim <- function(x=NULL, y=NULL, n.sims=2500, ns=seq(20, 200, 20), fake=TRUE){
    fr.rej <- 0
    k <- 0
    for (j in ns){
      k <- k+1
      fr.rej[k] <- 0
      for (i in 1:n.sims){
            if(fake){
              x <- runif(j)
              y <- rpois(j, 5)
            }
            fr.rej[k] <- fr.rej[k] +
              (perm(numperm = 2500, grad=sample(x, length(x)),
                  spAbun=y)[1]<0.05)/n.sims
        if (i%%100==0) print(paste("Sample size=", j,  ":", i, "of", n.sims, fr.rej[k]))
      }
    }
    return(fr.rej)
}

evergTP <- read.table(paste(dataDIR, "glades.env.txt", sep="/"), header=T)
evergInv <- read.table(paste(dataDIR, "glades.taxa.txt", sep="/"), header=T)
evergTP$Site<-row.names(evergTP)
names(evergTP)[1]<-"TP"
evergInv$Site <- row.names(evergInv)
glades.data <- merge(evergTP, evergInv, by="Site")

prI <- 0
j <- 0
for (i in names(glades.data)[-c(1,2)]){
    j <- j+1
    print(paste(i, j, "of", dim(glades.data)[2]-2, "taxa"))
    prI[j] <- typeIsim(x=glades.data$TP, y=glades.data[,i], n.sims=5000, ns=1, fake=F)
}
save(prI, file="PrTypeIGlades.RData")
prIsim <- typeIsim(n.sims=5000)
save(prIsim, file="PrTypeIN.RData")

load("PrTypeIGlades.RData")
load("PrTypeIN.RData")

tikz(file=paste(plotDIR, "typeIfig2.tex", sep="/"),
     height=3, width=4.5, standAlone=F)
par(mar=c(3,3,0.5,0.125), mgp=c(1.75,0.125,0), las=1, tck=0.01)
plot(prI ~ jitter(rep(dim(glades.data)[1], length(prI))),
       xlim=c(20,200), ylim=c(0.15,0.8),
     xlab="Sample Size", ylab="Pr(type I error)", col="gray")
lines(seq(20,200,20), prIsim)
dev.off()



zeros <- apply(glades.data[,-c(1,2)], 2, function(x) sum(x==0))
##pdf(file=paste(plotDIR, "zeroes.pdf", sep="/"), height=3, width=4.5)
tikz(file=paste(plotDIR, "titan0s.tex", sep="/"), height=3, width=4.5)
par(mar=c(3,3, 0.5,0.25), mgp=c(1.75,0.125,0), las=1, tck=0.01)
plot(zeros, prI, xlab="\\# of zeroes", ylab="Pr(type I error)", ylim=range(prIsim, prI))
dev.off()
