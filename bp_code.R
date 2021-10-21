#################################################################################
# r code for reproduction of MS1 graphs                                         #  
# MS1: Guidelines for Use of the Approximate Beta-Poisson Dose-Response Model   #
# April 2016           (Gang Xie (John) @ Smart Water Research Centre)                                                         #
#################################################################################
#------------------------------------------------
# For producing Fig. 2, Fig. 3, and Fig. 4, you need the following R package, and 
#  special R functions written by the primary author of the manuscript.

if("gsl" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("gsl"); require(gsl)} else{require(gsl)}  # you need this R library


# 
# Parameter estimation for Beta-Poisson Dose-Response models based on
# the deviance statistic defined by (7-35) (page 280, Haas1999).
#    New version
# 
newBPoptim <-  function(DD, NN, nn, inivalue=c(0.5,1),shapeA=c(0.0001,100), 
                        shapeB=c(0.01,500), showresult=TRUE, approxi=TRUE) {
  
  # This is a function for finding the optimal estimate of parameters
  # for a Beta-Poisson Dose-Response model  using 'optim' function in R
  # This  is for a 2-parameter specification distribution.
  # variable transformation may be needed for computation convergence reason.
  # 
  # DD: mean dose level, i.e. the actual number of microorganisms taken by a
  #     person Di ~ Poisson(DD).  
  # NN: number of people exposed to the microorganisms.
  # nn: number of people infected after exposure.
  # aa,bb: parameters of a Beta dist from which to generate 
  #     the rate or probability of infection by a single microorganism.
  
  nd = length(DD)
  # 
  aa <- log((inivalue[1] - shapeA[1])/(shapeA[2]-inivalue[1]))
  bb <- log((inivalue[2] - shapeB[1])/(shapeB[2]-inivalue[2])) 
  
  
  fr <- function(rx) {
    ry1 <- shapeA[1] + (shapeA[2]-shapeA[1])/(1+exp(-rx[1]))
    ry2 <- shapeB[1] + (shapeB[2]-shapeB[1])/(1+exp(-rx[2]))
    
    x1 <- ry1  # proxy for shapeA
    x2 <- ry2  # proxy for shapeB
    
    if(approxi==TRUE)      ppi = 1-(1+DD/x2)^(-x1)
    #
    if(approxi==FALSE)  {
      
      ppi  = 1- hyperg_1F1(x1,(x1+x2), -DD, give=FALSE, strict=TRUE)
      # A better way to calculate the hypergeometric functions, namely 
      #  1F1 and 2F1 is to use the R functions  'hyperg_1F1' and
      #  'hyperg_2F1' (need to be modified) in R package 'gsl'.
      
      # 
      
    } # end of 'if(approxi==FALSE)'
    
    minRR = -2*sum(nn*(log(ppi+0.0000000001)-log(nn/NN +0.0000000001)) + 
                     (NN-nn)* (log(1-ppi +0.0000000001)-log(1-nn/NN+0.0000000001))) 
    minRR
    
  }	
  
  out.obj <-  optim(c(aa,bb),fr)
  miniObj <- out.obj$value
  isconverge <- out.obj$convergence
  espar <- out.obj$par
  
  param1 <- shapeA[1] + (shapeA[2]-shapeA[1])/(1+exp(-espar[1]))
  param2 <- shapeB[1] + (shapeB[2]-shapeB[1])/(1+exp(-espar[2]))
  
  if(isconverge > 0 && showresult==TRUE) { 
    cat("Warning: the optimization process did not converge! \n") }
  if(isconverge ==0 && showresult==TRUE) { 
    cat("Good results: The optimization process converges! \n") }
  
  N50 = param2*(2^(1/param1) - 1)
  N50 = round(N50,3)
  
  
  
  if(showresult==TRUE)  cat("The estimated results are:\n\n")
  list(estshaA=param1,estshaB=param2, minimum=miniObj, N50=N50)
  
}       # end of the function

# Scenario: rate of infection of a single microorganism r ~ beta(a,b)
#---------

betaRsampA <- function(Ns, Ds, a, b) {
  # Ns: number of people exposed to the microorganisms, Ns could be a vector
  # Ds: mean dose level, i.e. the actual number of microorganisms taken by a
  #     person Di ~ Poisson(Ds).  Like the variable Ns, Ds could be a vector.
  # a,b: parameters of a Beta dist from which to generate 
  #     the rate or probability of infection by a single microorganism.
  #
  
  
  nx <- length(Ns)
  nd <- length(Ds)
  if(nx!=nd) stop("Ns and Ds must be the same length!")
  yy = as.numeric(nx)
  Nsum = cumsum(Ns); si = c(1,(Nsum[-nx]+1))
  
  n1 = sum(Ns)  # need to model the process based on each individual
  Ds1 = rep(Ds,Ns)
  yy1 = Di = as.numeric(n1)
  
  for (i in 1:n1) {
    Di[i] = rpois(1,Ds1[i])
    if(is.na(Di[i])==TRUE && Ds1[i] >1000000) Di[i] = Ds1[i]
    rr = rbeta(1,a,b)  # one prob for each individual
    pinfe1 <- 1- pbinom(0,Di[i],rr)
    yy1[i] = rbinom(1,1,pinfe1) }
  for (j in 1:nx)  yy[j] = sum(yy1[si[j]:Nsum[j]])
  yy            
}  # end of function


#-----------------------

ApprxiBPci = function(Ns, Ds, a, b, simuN, Dsim, appxi = TURE, pic=TRUE) {
  
  # Ns = number of people at health risk
  # Ds = mean dose administerted
  # a, b = paramters alpha and beta of an approximate B-P model
  # simuN : number of simulation runs
  # Dsim : plotting positions in terms of mean dose
  
  aaa = bbb = NULL 
  for(i in 1:simuN) {
    #   
    #   sampi = betaRsampA(Ns=Ns, Ds=Ds, a=aa[1],b=ba[1]) 
    sampi = betaRsampA(Ns, Ds, a,b)
    #   
    a0 = a/2; b0 = b/5
    aL = a/100; aU = ifelse(a<10, a*100, a*20)
    bL = b/1000; bU = ifelse(b<1000, b*100, b*20)
    
    out = newBPoptim(Ds,Ns,sampi, inivalue= c(a0,b0), shapeA=c(aL,aU),shapeB=c(bL,bU),
                     showresult=F, appxi)  
    #     
    aaa = c(aaa,out$estshaA)
    bbb = c(bbb, out$estshaB)
  }
  
  
  bpP = NULL; NN50 = NULL
  for(i in 1:length(aaa))  {
    if(appxi == TRUE)  midP = 1- (1+ Dsim/bbb[i])^(-aaa[i])
    if(appxi == FALSE) midP = 1- hyperg_1F1(aaa[i],(bbb[i]+aaa[i]),-Dsim)
    #     
    bpP = c(bpP,midP) 
    mid50 = bbb[i]*(2^(1/aaa[i])-1)
    NN50 = c(NN50, mid50)} 
  
  bpP.m = matrix(round(bpP,4),ncol=length(Dsim),byrow=T)
  
  poiE = 1- exp(-Dsim)
  py = y/Ns
  
  Pband = NULL
  for(j in 1:length(Dsim)) {
    
    midb = as.numeric(quantile(bpP.m[,j],c(0.025,0.5,0.975),na.rm=TRUE))
    Pband = c(Pband, midb) }
  Pband.m = matrix(Pband,ncol=length(Dsim))
  
  p25a = Pband.m[1,]; p50a = Pband.m[2,]; p975a= Pband.m[3,]
  #
  
  if(pic== TRUE) {
    #   op= par(mfrow=c(1,1), mar=c(2.5,3,1.5,1),mgp=c(1.5,0.5,0))
    
    xnames = c("0.01", "1", "100", expression(10^4),expression(10^6) , expression(10^8), 
               expression(10^10)) 
    xxx = c(-4.6, 0, 4.6, 9.21, 13.81, 18.42, 23.03)
    plot(log(Dsim),p50a,type="n",ylim=c(0,1), 
         xlab="Mean Dose", xaxt="n", ylab="Probability of Infection")
    #
    axis(1, at=xxx, labels= xnames)
    if(appxi == FALSE)  {lines(log(Dsim),p25a,lty=3)
      lines(log(Dsim),p975a,lty=3)
      lines(log(Dsim),p50a,lwd=2)   }
    
    if(appxi == TRUE)  {lines(log(Dsim),p25a,lty=3, col=2)
      lines(log(Dsim),p975a,lty=3, col=2)
      lines(log(Dsim),p50a, col=2)   }
    
    lines(log(Dsim),poiE,col=2, lwd=2,lty=2) 
    points(log(Ds), py)
    # par(op)   
  }
  #
  invisible(list(p25c=p25a,p50c=p50a, p975c=p975a, NN50=NN50,
                 alpha=aaa, beta=bbb))
}

#-------------------------------