############################################################################
############# Central weighting Matrix W / "kernels" matrix ################


##### ICM and Zheng tests central matrix

### Compute the kernel matrix given x and a bandwidth h

wmat <- function(x, h ,ker='normal',knorm="sd",remove=FALSE)
{
  #   ICM and smoothing matrix
  #   If no bandwidth is provided for the function, h=1 is used  (no smoothing)
  #   The principle diagonal is replaced with zeros if remove = TRUE.

  n<-dim(x)[1]

  mat <- apply(x,1,function(e) x - e*ones(n,1))

  # kernel smoothing

  wmat <-  kstand(mat / h,ker=ker,knorm=knorm)/h;

  # principle diagonal of the matrix replaced with zeros
  # "leave-one-out"

  if (remove==TRUE) wmat <-  wmat - diag(diag(wmat))

  wmat

}

### Standardized kernel smoothing function

kstand <- function(x,ker='normal',knorm="sd")
{

  # Densities such as the integral of the square of the density is one
  # if knorm is sq or such that the sd is one if knorm is sd.

  if (ker=='normal') # Normal pdf
  {
    s <- switch(knorm,
                sd = 1,
                sq = 1/(2*sqrt(pi)))
    aux <- exp(- x^2 /(2*s^2)) / (sqrt(2*pi)*s)
  }

  if (ker=='triangle') # Triangle pdf
  {
    s <- switch(knorm,
                sd = sqrt(6),
                sq = (2/3))
    aux <- (s-abs(x))
    aux <- aux * (aux >0)/s^2
  }

  if (ker=='logistic') # Logistic pdf
  {
    s <- switch(knorm,
                sd = 1/sqrt(2),
                sq = 1/4)
    aux <- exp(-abs(x)/s)/(2*s)
  }

  if (ker=="sinc") # Sine-cardinal function
  {
    aux <- sqrt(3/2)*sinc(x)^2
    aux <- apply(aux, 1, function(x) replace(x,list=which(x=="NaN"),values=1))
  }

  aux

}

### Standardized central matrix

Wemult <- function(Z,h=1,ker="normal",remove=FALSE,knorm="sd")
{
  # Compute standardized Kernel matrix for multivariate Z

  k <- dim(Z)[2]
  n <- dim(Z)[1]

  if (k > 1)
  {
    if (sd(Z[,1])!=0){
      Ztnorm <- as.matrix(Z[,1]/sd(Z[,1]))
      We <- wmat(Ztnorm,h=h,ker =ker ,remove=remove,knorm=knorm)
    } else {
      We<-matrix(kstand(x=0,ker=ker,knorm=knorm),n,n)
      # In case there is an intercept
    }

    for (j in 2:k)
    {
      if (sd(Z[,j])!=0){
        Ztnorm <- as.matrix(Z[,j]/sd(Z[,j]))
      } else {
        Ztnorm<-as.matrix(Z[,j])
      }
      # the multivariate kernel function is just the product
      # of univariate kernel function
      We <- We * wmat(Ztnorm,h=h,ker=ker,remove=remove,knorm=knorm)
    }
  } else {
    We <- wmat(Z/sd(Z),h=h,ker=ker,remove=remove,knorm=knorm)
  }
  We
}


##### Escanciano test central matrix

# Similar to the ICM weighting matrix except that
# kernel functions are not used but different functions
# which simplify integration over the 1-hypersphere
# See Escanciano's paper for more details

### Individual elements of the matrix

WescI<-function(i,j,xr,x){

  k<-length(xr)
  if (k>1){
    if (x[i,]==xr || x[j,]==xr) {
      Aijr<-pi
    } else {
      Aij0<-abs(pi-acos(t(x[i,]-xr)%*%(x[j,]-xr)/(norm(x[i,]-xr,type="F")*norm(x[j,]-xr,type="O"))))
      Aijr<-Aij0*pi^(k/2-1)/gamma(k/2+1)
    }
  } else if (k==1){
    if (x[i]==xr || x[j]==xr) {
      Aijr<-pi
    } else {
      Aij0<-abs(pi-acos(t(x[i]-xr)*(x[j]-xr)/(norm(x[i]-xr,type="F")*norm(x[j]-xr,type="O"))))
      Aijr<-Aij0*pi^(k/2-1)/gamma(k/2+1)
    }
  }

  Aijr

}

### Creation of the matrix from half of its elements
### (the matrix is symmetric)

Wesc<-function(x,remove=F){
  n=dim(x)[1]

  WescI2<-function(i,j,x){
    Aij<-apply(x,1, function(e) WescI(i,j,as.matrix(e),x))
    sum(Aij)
  }

  WescI3<-function(i,x){

    n=dim(x)[1]
    sapply(seq(i,n), function(j) WescI2(i,j,x))
  }


  WW<-sapply(seq(1,n), function(i) WescI3(i,x))

  W<-matrix(0,n,n)

  for (i in 1:(n)){
    W[i,i:n]<-WW[[i]]
  }
  Wf<-W+t(W)


  if (remove==T){
    diag(Wf)<-0
  } else if (remove==F){
    diag(Wf)<-diag(Wf)/2
  }

  Wf

}


##### SICM test central matrix

# Here instead of integrating over the whole hypersphere
# sums are used and a direction ca be given

### Finds a random beta in the hypersphere norm 1
### which may follow the direction given initially

Randhyper<-function(direct){

  Ok<-F

  while (Ok==F){

    b1<-as.matrix(rnorm(length(direct)))
    b2<-norm(b1,type="F")
    b3<-b1/b2

    s<-sum(sign(b3))
    sd<-sum(direct)
    sd0<-sum(direct==0)

    if ( s <= sd+sd0 & s >=sd-sd0){

      Ok=T
      bf<-sortr(b3,direct)


    } else if ( -s<= sd+sd0 & -s>=sd-sd0){

      Ok=T
      bf<-sortr(-b3,direct)

    }
  }

  bf

}

### Reorders the random beta so that it follows
### the initial direction

sortr<-function(x,direct){

  d1<-which(direct==1)
  ld1<-length(d1)
  dm1<-which(direct==-1)
  ldm1<-length(dm1)

  x1<-which(x>0)
  lx1<-length(x1)
  xm1<-which(x<0)
  lxm1<-length(xm1)

  b0<-rep(0,length(direct))

  if (length(dm1)==0 & length(d1)==0){
    b0<-x
  } else if (length(dm1)==0){
    b0[d1]<-x[x1[1:ld1]]
    b0[b0==0]<-x[-x1[1:ld1]]
  } else if (length(d1)==0){
    b0[dm1]<-x[xm1[1:ldm1]]
    b0[b0==0]<-x[-xm1[1:ldm1]]
  } else {
    b0[dm1]<-x[xm1[1:ldm1]]
    b0[d1]<-x[x1[1:ld1]]
    b0[b0==0]<-x[-c(xm1[1:ldm1],x1[1:ld1])]
  }

  b0

}

### Builds the individual elements of the weighting matrix

Ws1<-function(i,j,x,hyper,h=1,ker="normal",knorm="sd"){

  n<-dim(x)[1]
  k<-dim(x)[2]

  if (i==j){
    Wij<-dim(hyper)[2]*kstand(0,ker=ker,knorm=knorm)
  } else if (i!=j){
    if (k==1){
      Wijb<-sapply(hyper, function(e) kstand((x[i]-x[j])*e/h,ker=ker,knorm=knorm))
    } else if (k>1){
      if (sd(x[,1])!=0){
        Wijb<-apply(hyper,2, function(e) kstand((x[i,]-x[j,])%*%e/h,ker=ker,knorm=knorm))
      } else if (sd(x[,1])==0){
        Wijb<-apply(hyper,2, function(e) kstand((x[i,2:k]-x[j,2:k])%*%e/h,ker=ker,knorm=knorm))
      }
    }

    Wij<-sum(Wijb)

  }

  Wij

}

### Builds the matrix from helf its elements

Ws3<-function(x,nbeta=200,h=1,ker="normal",knorm="sd",remove=F,direct){

  n<-dim(x)[1]
  k<-dim(x)[2]

  # Preliminary standardization

  if (sd(x[,1])==0){
    x<-cbind(1,apply(x[,2:k],2,function(e) e/sd(e)))
  } else {
    x<-apply(x,2,function(e) e/sd(e))
  }

  hyper<-sapply(1:nbeta, function(e) Randhyper(direct))


  Ws2<-function(i,x,hyper,h,ker,knorm){

    n<-dim(x)[1]

    W1<-sapply(i:n, function(e) Ws1(i=i,e,x=x,hyper=hyper,h=h,ker=ker,knorm=knorm))
    W1

  }


  W2<-sapply(1:n, function(e) Ws2(e,x=x,hyper=hyper,h=h,ker=ker,knorm=knorm))

  Wf<-matrix(0,n,n)
  for (i in 1:n){
    Wf[i,i:n]<-W2[[i]]
  }
  Wf<-Wf+t(Wf)

  if (remove==T){
    diag(Wf)<-0
  } else if (remove==F){
    diag(Wf)<-diag(Wf)/2
  }

  Wf/nbeta

}


##### PALA test weighting matrix

# Similar to SICM except that the statistic
# uses the random beta which maximizes it

### Finds a random beta in the hypersphere of norm 1

Randsimple<-function(k){

  b1<-as.matrix(rnorm(k))
  b2<-norm(b1,type="F")
  b3<-b1/b2
  b3

}

### Build a kernel matrix for a given beta

Wpalasingle<-function(x,beta,h=1,ker='normal',knorm="sd",remove=T){

  k<-dim(x)[2]
  n<-dim(x)[1]

  if (sd(x[,1])==0){
    x<-cbind(1,apply(x[,2:k],2,function(e) e/sd(e)))
  } else {
    x<-apply(x,2,function(e) e/sd(e))
  }

  if (k>1){
    xbeta<-x%*%as.matrix(beta)
  } else if (k==1){
    xbeta<-x*beta
  }

  Wbeta<-wmat(x=xbeta,h=h,ker=ker,knorm=knorm,remove=remove)
  Wbeta

}

### The matrices of beta are then stuck together

Wpalafull<-function(x,h=1,beta0,nbeta=200,ker='normal',knorm="sd",remove=T){

  n<-dim(x)[1]
  k<-dim(x)[2]

  hyper<-sapply(1:nbeta, function(e) Randsimple(k))

  if (k>1){
    Wf<-apply(hyper,2, function(e) Wpalasingle(x=x,e, h=h ,ker=ker,knorm=knorm,remove=remove))
    Wb<-apply(hyper,2, function(e) prod(e!=beta0))
  } else if (k==1){
    Wf<-sapply(hyper, function(e) Wpalasingle(x=x,e, h=h ,ker=ker,knorm=knorm,remove=remove))
    Wb<-sapply(hyper, function(e) prod(e!=beta0))
  }

  rbind(Wf,Wb,rep(h,nbeta))

}



#########################################################################
########### Non-parametric estim of variance of errors ################

# Yin, Geng, Li and Wang (2010) estimator
# Y is typically the vector of residuals
# We is the matrix of kernel built from the regressors

Vest <- function(Y,We)
{
  Wsum <- apply(We,1,sum)
  f1 <- (We%*%Y)/Wsum
  (We%*%((Y-f1)^2))/Wsum
}

