################################################################
#################### Individual Statistics #####################

### ICM, SICM and ESCA statistic

mystat <- function(uhat,mat)
{
  n <- length(uhat)
  s <- (t(uhat)%*%mat%*%uhat)/n
  as.numeric(s)
}

### Zheng statistic

myzheng <- function(uhat,mat,matt)
{
  n <- length(uhat)
  s <- (t(uhat)%*%mat%*%uhat)/(n*(n-1))
  var_cond <- Vest(uhat,matt)
  mat2 <- mat^2
  v2 <- 2*t(var_cond)%*%mat2%*%var_cond/(n^2*(n-1)^2)
  as.numeric(s /sqrt(v2))
}

### PALA statistic

mypala <- function(uhat,mat,matt)
{
  n <- length(uhat)
  nbeta<-dim(mat)[2]
  var_cond <- Vest(uhat,matt)
  s<-sapply(1:nbeta, function(e) sqrt(mat[(n^2+2),1])*(uhat)%*%matrix(mat[1:n^2,e],n,n)%*%uhat/(n-1)-log(n)*sqrt(2*t(var_cond)%*%matrix(mat[1:n^2,e],n,n)%*%var_cond*mat[(n^2+1),e]/((n-1)*sqrt(mat[(n^2+2),1])))/n^(1.5))
  as.numeric(max(s))

}



########################################################################
###################### Test statistic and  pvalue ######################


##### Function which guards against some user's potential mistakes
##### Not that important

SpeTest<-function(eq="default",x="default",y="default",inter=TRUE,type="icm",boot="wild",nboot=99,ker="normal",knorm="sd",hv="default",cch="default",nbeta="default",beta0="default",direct="default")
{

  if (class(x)!="character" & class(y)!="character" & class(eq)=="character"){
    if (inter==T){
      x<-cbind(1,x)
    }

    StatPval(x=x,y=y,nboot=nboot,hv=hv,cch=cch,nbeta=nbeta,beta0=beta0,ker=ker,knorm=knorm,boot=boot,type=type,direct=direct)

  } else if (class(x)=="character" & class(y)=="character" & class(eq)!="character"){

    if ((names(eq$coefficient[1])=="(Intercept)")==T){
      x<-cbind(1,eq$model[,2])
    } else if ((names(eq$coefficient[1])=="(Intercept)")==F) {
      x<-eq$model[,2]
    }

    y<-eq$model[,1]

    StatPval(x=x,y=y,nboot=nboot,hv=hv,cch=cch,nbeta=nbeta,beta0=beta0,ker=ker,knorm=knorm,boot=boot,type=type,direct=direct)

  } else {
    cat(" \n  Provide either a model of class lm or a dependent variable and regressors. \n ")
  }

}


##### Main function which delivers the statistic and its pvalue

StatPval<-function(x,y,type="icm",boot="wild",nboot=99,ker="normal",knorm="sd",hv="default",cch="default",nbeta=99,beta0="default",direct="default")
{

  ##### To begin with

  x<-as.matrix(x)
  eq<-lm(y~x+0)
  uhat <- residuals(eq)
  fit <- fitted(eq)

  n<-length(fit)
  k<-dim(x)[2]

  ##### Options

  ### "kernels" matrix bandwidth for Zheng, PALA and SICM
  ### By default,

  if (cch=="default"){
    if (type=="zheng"){
      cch=1.06*n^(-1/(4+k))
    } else if (type=="sicm" || type=="pala"){
      cch=1.06*n^(-1/5)
    }
  }

  ### Non-parametric conditional variance estimator bandwidth
  ### Default is rule of thumb scaled by parameter dimension

  if (hv=="default") {
    hv=1.06*n^(-1/(4+k))
  }

  ### "kernels" matrix leave one out or not
  ### default depends on each statistic general

  if (type=="zheng"||type=="sicm"||type=="pala"){
    remove=T
  } else if (type=="icm"||type=="esca"){
    remove=F
  }

  ### In SICM, initial direction for the parameter can be given

  if (class(direct)=="character"){
    if (sd(x[,1])==0){
      direct<-rep(0,k-1)
    } else {
      direct<-rep(0,k)
    }
  }

  ### In PALA initial prediction for the parameter can be given

  if (beta0=="default"){
    beta0<-eq$coefficients
  }

  ### In PALA and SICM by default the number of beta in the
  ### hypersphere increases with the regressors dimension

  if (nbeta=="default"){
    nbeta=200*sqrt(k)
  }

  ##### Statistics' central matrix which depends on the type

  W<-eye(n)
  We<-W

  if (type=="icm"){
    W<-Wemult(x,h=1,ker=ker,remove=remove,knorm=knorm)
    if (boot=="smooth"){
      We<-Wemult(x,h=hv,ker=ker,remove=remove,knorm=knorm)
    }
  } else if (type=="zheng"){
    W<-Wemult(x,h=cch,ker=ker,remove=remove,knorm=knorm)
    if (cch==hv){
      We<-W
    } else {
      We<-Wemult(x,h=hv,ker=ker,remove=remove,knorm=knorm)
    }
  } else if (type=="esca"){
    W<-Wesc(x,remove=remove)
    if (boot=="smooth"){
      We<-Wemult(x,h=hv,ker=ker,remove=remove,knorm=knorm)
    }
  } else if (type=="sicm"){
    W<-Ws3(x,nbeta=nbeta,h=cch,ker=ker,remove=remove,knorm=knorm,direct=direct)
    if (boot=="smooth"){
      We<-Wemult(x,h=hv,ker=ker,remove=remove,knorm=knorm)
    }
  } else if (type=="pala"){
    W<-Wpalafull(x,h=cch,beta0=beta0,nbeta=nbeta,ker=ker,knorm=knorm,remove=remove)
    We<-Wemult(x,h=hv,ker=ker,remove=remove,knorm=knorm)
  }


  ##### The statistic is computed

  SStat <- switch(type,"icm"=mystat(uhat=uhat,mat=W),"sicm"=mystat(uhat=uhat,mat=W),"esca"=mystat(uhat=uhat,mat=W),"zheng"=myzheng(uhat=uhat,mat=W,matt=We),"pala"=mypala(uhat=uhat,mat=W,matt=We))

  ##### Bootstraps

  Sboot<-fboot(X=x,uhat=uhat,fit=fit,mat=W,matt=We,nboot=nboot,boot=boot,type=type)

  ##### P-values

  PPval <- mean(Sboot>SStat)

  ##### Recap

  x<-list(SStat,PPval,type,boot,nboot,ker,knorm,hv,cch,nbeta,beta0,direct)
  names(x)<-c("stat","pval","type","boot","nboot","ker","knorm","hv","cch","nbeta","beta0","direct")

  class(x)<-"STNP"
  return(x)

}



###########################################################################
########################### Presentation functions ########################

##### Print function

print.STNP<-function(x){

  name<-switch(x$type,"zheng"="Zheng test","icm"="Integrated Conditional Moment test","sicm"="Smooth Integrated Conditional Moment test","esca"="Escanciano test","pala"="Lavergne and Patilea test")

  if (x$type=="zheng"){
    pvalasy <- pnorm(x$stat,lower.tail = FALSE)
    statf <- c("Statistic : ",round(x$stat,digits=5))
    apvalf<-c("Asymptotic p-value : ",round(pvalasy,digits=5))
    pvalf<-c("Bootstrap p-value : ",round(x$pval,digits=5))

    cat("\n ",name,"\n\n ",statf,"\n ",apvalf,"\n ",pvalf,"\n ")

  } else {
    statf <- c("Statistic : ",round(x$stat,digits=5))
    pvalf<-c("Bootstrap p-value : ",round(x$pval,digits=5))

    cat("\n ",name,"\n\n ",statf,"\n ",pvalf,"\n ")

  }

}


##### Summary function

summary.STNP<-function(x){

  name<-switch(x$type,"zheng"="Zheng test","icm"="Integrated Conditional Moment test","sicm"="Smooth Integrated Conditional Moment test","esca"="Escanciano test","pala"="Lavergne and Patilea test")

  if (x$type=="zheng"){
    pvalasy <- pnorm(x$stat,lower.tail = FALSE)
    statf <- c("Statistic : ",round(x$stat,digits=5))
    apvalf<-c("Asymptotic p-value : ",round(pvalasy,digits=5))
    pvalf<-c("Bootstrap p-value : ",round(x$pval,digits=5))

    cat("\n ",name,"\n\n ",statf,"\n ",apvalf,"\n ",pvalf)

  } else {
    statf <- c("Statistic : ",round(x$stat,digits=5))
    pvalf<-c("Bootstrap p-value : ",round(x$pval,digits=5))

    cat("\n ",name,"\n\n ",statf,"\n ",pvalf)

  }

  bootp<-switch(x$boot,"wild"="Bootstrap type :  Wild bootstrap","smooth"="Bootstrap type :  Smooth Conditional Moments bootstrap")

  nbootp<-c("Number of bootstrap samples : ",x$nboot)
  kerp<-switch(x$ker,"normal"="Kernel function :  Normal p.d.f","triangle"="Kernel function :  Triangular p.d.f","logistic"="Kernel function :  Logistic p.d.f","sinc"="Kernel function :  Sine Cardinal function")
  knormp<-switch(x$knorm,"sd"="Standardized Kernel :  Standard deviation using Kernel density = 1","sq"="Standardized Kernel :  Integral of squared Kernel function = 1")

  cat("\n\n ",bootp,"\n ",nbootp,"\n\n ",kerp,"\n ",knormp)

  if (x$type=="pala"){
    beta0p<-c("Initial prediction for beta = ",x$beta0)
    cat("\n\n ",beta0p)
  }

  if (x$type=="sicm" & sum(x$direct!=0)>0){
    directp<-c("Initial direction for beta = ",round(x$direct,digits=5))
    cat("\n\n ",directp)
  }

  if (x$boot=="smooth" || x$type=="zheng" ||x$type=="pala"){
    hvp<-c("Non-parametric conditional variance estimator bandwidth = ",round(x$hv,digits=5))
    cat("\n ",hvp)
  }

  if (x$type!= "esca" & x$type!= "icm"){
    cchp<-c("Non-parametric statistic central kernel matrix bandwidth = ",round(x$cch,digits=5))
    cat("\n ",cchp)
  }

  cat(" \n ")

}

