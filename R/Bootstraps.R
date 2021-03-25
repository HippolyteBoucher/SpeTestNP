#############################################################################
############################# Bootstraps ####################################


### Recovers bootstrap statistic using initial regression

# the matrix w is random and is a faster to create the
# bootstraps' statistics

wboot <- function(w,e,mat,matt,us,fit,type)
{

  y <- fit + w*us
  eq <- lm(y~e)
  uhatb <- residuals(eq)
  switch(type,"icm"=mystat(uhatb,mat),"sicm"=mystat(uhatb,mat),"esca"=mystat(uhatb,mat),"zheng"=myzheng(uhatb,mat,matt),"pala"=mypala(uhatb,mat,matt))

}

### Builds the bootstrap statistics vector

fboot <- function(X,uhat,fit,mat,matt,nboot,boot,type)
{

  # w is created and standardized

  n <- length(uhat)
  pr <- 0.72360679774998   # (5+sqrt(5))/10
  c1 <-  - 0.6180339887499   # (1-sqrt(5))/2
  c2 <-  1.6180339887499    # (1+sqrt(5))/2
  vb <- runif(n*nboot)
  d1 <- (vb < pr)
  w <- matrix((c1*d1)+(c2*(1-d1)),nrow=n)

  ### Wild bootstraps or Smooth Conditional Moments bootstraps

  # What is redrawn is either the residual itself or the estimated
  # standard deviation

  u<-switch(boot,"wild"=uhat,"smooth"=sqrt(Vest(uhat,matt)))

  b <- apply(w,2, function(l) wboot(l,e=X,mat=mat,matt=matt,us=u,fit=fit,type=type))
  b

}

