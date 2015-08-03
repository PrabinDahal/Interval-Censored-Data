###########################################
# Dealing with Interval-Censored Data
###########################################
rm(list=ls())
library(survival);library(KMsurv);library(MIICD)

data(bcdeter)
dat<-bcdeter

# Convert right censored limit to infinity
dat$upper[is.na(dat$upper)]<-Inf
dat$cens<-ifelse(is.finite(dat$upper), 1, 0)
colnames(dat)<-c("left","right","treat","cens")
# Only radiotherapy arm kept for analysis
dat<-dat[which(dat$treat==1),] 

#---------------------------------------------------------------
# TURNBULL's ALGORITHM (Written by: Prof. Suely Ruiz Giolo)
# Available from: http://www.est.ufpr.br/rt/Turnbull.R
# Access Date: 03-08-2015
#---------------------------------------------------------------
cria.tau <- function(data){
  l <- data$left
  r <- data$right
  tau <- sort(unique(c(l,r[is.finite(r)])))
  return(tau)
}

S.ini <- function(tau){
  m<-length(tau)
  ekm<-survfit(Surv(tau[1:m-1],rep(1,m-1))~1)
  So<-c(1,ekm$surv)
  p <- -diff(So)
  return(p)
}

cria.A <- function(data,tau){
  tau12 <- cbind(tau[-length(tau)],tau[-1])
  interv <- function(x,inf,sup) ifelse(x[1]>=inf & x[2]<=sup,1,0)
  A <- apply(tau12,1,interv,inf=data$left,sup=data$right)
  id.lin.zero <- which(apply(A==0, 1, all))
  if(length(id.lin.zero)>0) A <- A[-id.lin.zero, ] 
  return(A)
}

Turnbull <- function(p, A, data, eps=1e-5,
                     iter.max=1000, verbose=TRUE){
  n<-nrow(A)
  m<-ncol(A)
  Q<-matrix(1,m)
  iter <- 0
  repeat {
    iter <- iter + 1
    diff<- (Q-p)
    maxdiff<-max(abs(as.vector(diff)))
    if (verbose)
      print(maxdiff)
    if (maxdiff<eps | iter>=iter.max)
      break
    Q<-p
    C<-A%*%p
    p<-p*((t(A)%*%(1/C))/n)
  }
    cat("Iterations = ", iter,"\n")
    cat("Max difference = ", maxdiff,"\n")
    cat("Convergence criteria: Max difference < 1e-5","\n")
  dimnames(p)<-list(NULL,c("P Estimate"))
  surv<-round(c(1,1-cumsum(p)),digits=7)
  right <- data$right
   if(any(!(is.finite(right)))){
    t <- max(right[is.finite(right)])
    return(list(time=tau[tau<t],surv=surv[tau<t]))
  }
  else
    return(list(time=tau,surv=surv))
}

########################################
# Comparison of three different methods
########################################

#----------------------------------
# Method I: Ignorant Surv Estimates
#----------------------------------
dat$naive<-ifelse(is.finite(dat$right), dat$right, dat$left)
naive.surv<-survfit(Surv(naive,cens)~1, data=dat)

#------------------------------------------------
# Method II: Turnbull's Interval Censored Method
#------------------------------------------------
tau<-cria.tau (dat)
p<-S.ini(tau)
A<-cria.A(dat,tau)
turn.surv<-Turnbull (p, A, dat)
int.cens<-data.frame(time=turn.surv[1], surv=turn.surv[2])

#---------------------------------------------------------------
# Method III: Data Augmentation and Multiple Imputation Approach
#----------------------------------------------------------------
set.seed(1)
res1<-DA.surv(k = 1000, m = 50, data = dat, conf.int = TRUE , alpha = 0.05 )
int<-res1$est

#-------------------------------------
# Plot the three estimates together
#-------------------------------------
plot(naive.surv,  # Method I
	xlab="Months on study", 
	ylab="Estimated Survival Function",
	conf.int=FALSE, las=1, lwd=2, 
	cex.lab=1.5, cex.axis=1.2, col="royalblue3"
)
lines(int$surv~int$time, type ="s", ylim=c(0.9,1), col="springgreen4",lwd=2) # Method II
lines(int.cens$surv~int.cens$time, type="s", col="red", lwd=2) # Method-III

legend("bottomleft",c("Turnbull's Method for Interval Censored Data","Data Augmentation and Multiple Imputation Approach","Ignorant Kaplan-Meier"), 
	col=c("red","springgreen4","royalblue3"),
	lty=1, bty="n", lwd=2)


# Session Infro
sessionInfo()


