#### Poisson likelihood and NR scheme ####

#### 2.7 Poisson GLM: using fit701 function ####
llmaxM2B=function(b1,b2,b3,k2,g3,dv,ev,wv=1){
  #   b1,b3,k2,g3 are given
  #   solve for b2
  b21=b2
  b20=b21-1
  thetat=k2*ev*exp(b1+b3*g3)
  s1=sum(dv*k2*wv)
  while(abs(b21-b20) > 0.1)
  {
    b20=b21
    f0=sum((exp(b20*k2)*thetat)*wv)-s1
    df0=sum((exp(b20*k2)*k2*thetat)*wv)
    b21=b20-f0/df0
  }
  b21
}
llmaxM2C=function(b1,b2,b3,k2,g3,dv,ev,wv=1){
  #   b1,b2,k2,g3 are given
  #   solve for b3
  b31=b3
  b30=b31-1
  thetat=g3*ev*exp(b1+b2*k2)
  s1=sum(dv*g3*wv)
  while(abs(b31-b30) > 0.1)
  {
    b30=b31
    f0=sum((exp(b30*g3)*thetat)*wv)-s1
    df0=sum((exp(b30*g3)*g3*thetat)*wv)
    b31=b30-f0/df0
  }
  b31
}


llmaxM2D=function(b1,b2,b3,k2,g3,dv,ev,wv=1){
  #   b1,b2,b3,g3 are given
  #   solve for k2
  k21=k2
  k20=k21-1
  thetat=b2*ev*exp(b1+b3*g3)
  s1=sum(dv*b2*wv)
  while(abs(k21-k20) > 0.1)
  {
    k20=k21
    f0=sum((exp(k20*b2)*thetat)*wv)-s1
    df0=sum((exp(k20*b2)*b2*thetat)*wv)
    k21=k20-f0/df0
  }
  k21
}



llmaxM2E=function(b1,b2,b3,k2,g3,dv,ev,wv=1){
  #   b1,b2,b3,k2 are given
  #   solve for g3
  g31=g3
  g30=g31-1
  thetat=b3*ev*exp(b1+b2*k2)
  s1=sum(dv*b3*wv)
  while(abs(g31-g30) > 0.1)
  {
    g30=g31
    f0=sum((exp(g30*b3)*thetat)*wv)-s1
    df0=sum((exp(g30*b3)*b3*thetat)*wv)
    g31=g30-f0/df0
  }
  g31
}


logit=function(x){ log(x/(1-x)) }


llmaxM7A=function(k1,k2,k4,b2,b3,b4,g3,dv,ev,wv=1){
  # Optimise over kappa1(t) given kappa2(t), kappa4(t), beta2(x), beta3(x), beta4(x) and gamma3(t-x)
  h=0.000001
  mv=log(1+exp(k1+k2*b2+k4*b4+g3*b3))
  dm=1/(1+exp(-k1-k2*b2-k4*b4-g3*b3))
  f0=sum((dv/mv*dm-ev*dm)*wv)
  mvh=log(1+exp((k1+h)+k2*b2+k4*b4+g3*b3))
  dmh=1/(1+exp(-(k1+h)-k2*b2-k4*b4-g3*b3))
  f0h=sum((dv/mvh*dmh-ev*dmh)*wv)
  df0=(f0h-f0)/h
  k1new=k1-f0/df0
  k1new
}

llmaxM7B=function(k1,k2,k4,b2,b3,b4,g3,dv,ev,wv=1){
  # Optimise over kappa2(t) given kappa1(t), kappa4(t), beta2(x), beta3(x), beta4(x) and gamma3(t-x)
  h=0.000001
  mv=log(1+exp(k1+k2*b2+k4*b4+g3*b3))
  dm=b2/(1+exp(-k1-k2*b2-k4*b4-g3*b3))
  f0=sum((dv/mv*dm-ev*dm)*wv)
  mvh=log(1+exp(k1+(k2+h)*b2+k4*b4+g3*b3))
  dmh=b2/(1+exp(-k1-(k2+h)*b2-k4*b4-g3*b3))
  f0h=sum((dv/mvh*dmh-ev*dmh)*wv)
  df0=(f0h-f0)/h
  k2new=k2-f0/df0
  k2new
}


llmaxM7D=function(k1,k2,k4,b2,b3,b4,g3,dv,ev,wv=1){
  # Optimise over kappa4(t) given kappa1(t), kappa2(t), beta2(x), beta3(x), beta4(x) and gamma3(t-x)
  h=0.000001
  mv=log(1+exp(k1+k2*b2+k4*b4+g3*b3))
  dm=b4/(1+exp(-k1-k2*b2-k4*b4-g3*b3))
  f0=sum((dv/mv*dm-ev*dm)*wv)
  mvh=log(1+exp(k1+k2*b2+(k4+h)*b4+g3*b3))
  dmh=b4/(1+exp(-k1-k2*b2-(k4+h)*b4-g3*b3))
  f0h=sum((dv/mvh*dmh-ev*dmh)*wv)
  df0=(f0h-f0)/h
  k4new=k4-f0/df0
  k4new
}


llmaxM7C=function(k1,k2,k4,b2,b3,b4,g3,dv,ev,wv=1){
  # Optimise over gamma3(t-x) given kappa1(t), kappa2(t), kappa4(t), beta2(x), beta3(x), beta4(x) 
  h=0.000001
  mv=log(1+exp(k1+k2*b2+k4*b4+g3*b3))
  dm=b3/(1+exp(-k1-k2*b2-k4*b4-g3*b3))
  f0=sum(wv*(dv/mv*dm-ev*dm))
  mvh=log(1+exp(k1+k2*b2+k4*b4+(g3+h)*b3))
  dmh=b3/(1+exp(-k1-k2*b2-k4*b4-(g3+h)*b3))
  f0h=sum(wv*(dv/mvh*dmh-ev*dmh))
  df0=(f0h-f0)/h
  g3new=g3-f0/df0
  g3new
}


llmaxM7F=function(x0,k1,k2,k4,b2,b4,g3,dtx,etx,xv,wa=1){
  #   b2,b4,k1,k2,k4,g3 are given
  #   solve for x0 which dictates the form of beta3
  x01=x0
  x00=x01-1000
  h=0.000001
  m=length(k1)
  n=length(b2)
  k1a=array(k1,c(m,n))
  k2a=array(k2,c(m,n))
  k4a=array(k4,c(m,n))
  b2a=t(array(b2,c(n,m)))
  b4a=t(array(b4,c(n,m)))
  b3=x0-xv
  b3a=t(array(b3,c(n,m)))
  g3a=b3a*0
  for(i in 1:m)
  {
    g3a[i,]=g3[(n+i-1):i]
  }
  while(abs(x01-x00) > 100)
  {
    x00=x01
    mv=log(1+exp(k1a+k2a*b2a+k4a*b4a+g3a*b3a))
    dm=g3a/(1+exp(-k1a-k2a*b2a-k4a*b4a-g3a*b3a))
    f0=sum(wa*(dtx/mv*dm-etx*dm))
    mvh=log(1+exp(k1a+k2a*b2a+k4a*b4a+g3a*(b3a+h)))
    dmh=g3a/(1+exp(-k1a-k2a*b2a-k4a*b4a-g3a*(b3a+h)))
    f0h=sum(wa*(dtx/mvh*dmh-etx*dmh))
    df0=(f0h-f0)/h
    x01=x00-f0/df0
  }
  x01
}
fit701=function(xv,yv,etx,dtx,wa){
  # Model M1
  # Lee-Carter model
  # log m(t,x) = beta1(x) + beta2(x).kappa2(t) + Poisson error
  
  # Inputs:
  #   xv = vector of ages, length n
  #   yv = vector of years, length m
  #   etx = m x n matrix of exposures
  #   dtx = m x n matrix of deaths
  #   wa = m x n matrix of weights (0 or 1)
  xv<-as.vector(unlist(xv))
  yv<-as.vector(unlist(yv))
  etx<-as.matrix(etx)
  dtx<-as.matrix(dtx)
  wa<-as.matrix(wa)
  mtx=dtx/etx	  # matrix of death rates
  
  qtx=1-exp(-mtx) # matrix of mortality rates
  
  #if(max(xv) > 89)
  #{  
  #  cat("Upper age too high - suggest abort programme\n") 
  #} # not applicable for our dataset, revision by Sun
  
  n=length(xv)	# number of ages
  m=length(yv)	# number of years
  
  cy=(yv[1]-xv[n]):(yv[m]-xv[1])  # cohort approximate years of birth
  
  # initialise parameter vectors
  beta1v=(1:n)*0
  beta2v=(1:n)*0
  beta3v=(1:n)*0		# dummy vector, this will stay at 0
  kappa2v=(1:m)*0
  gamma3v=(1:(n+m-1))*0	# dummy vector, this will stay at 0
  ia=array((1:m),c(m,n))	# matrix of year indexes, i, for the data
  ja=t(array((1:n),c(n,m)))	# matrix of age indexes, j, for the data
  ya=ia-ja		 	# matrix of year of birth indexes for the data
  imj=(1-n):(m-1)		# the range of values taken by i-j
  lg=n+m-1		 	# number of different values taken by i-j
  ca=ya+yv[1]-xv[1]		# matrix of years of birth
  
  # Now set weights to zero for cohorts with fewer than 5 observations
  for(k in 1:lg)
  {
    nk=sum((ca == cy[k])*wa)
    if(nk < 5)
    {
      wa=wa*(1- (ca == cy[k]))
    }
  }
  
  ww=cy*0+1	 # this is a vector of 1's and 0's with
  # a 0 if the cohort is completely excluded
  for(k in 1:lg)
  {
    ww[k]=ww[k]*(sum((ca == cy[k])*wa) > 0)
  }
  
  # Stage 0
  # Gives initial estimates for beta1(x), beta2(x) and kappa2(t)
  mx=mean(xv)
  for(j in 1:n)
  {
    beta1v[j]=sum(log(mtx[,j])*wa[,j])/sum(wa[,j])
    beta2v[j]=1/n
  }
  kappa2v=(m:1)-(m+1)/2
  
  # Stage 1: iterate
  l0=-1000000
  l1=-999999
  iteration=0
  # l1 is the latest estimate of the log-likelihood
  # l0 is the previous estimate
  # we continue to iterate if the improvement in log-likelihood
  # exceeds 0.0001
  while(abs(l1-l0) > 0.0001)
  {
    iteration=iteration+1
    
    l0=l1
    # Stage 1B optimise over the beta2(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      dv=dtx[,j]	# actual deaths
      ev=etx[,j]	# exposure
      beta2v[j]=llmaxM2B(beta1v[j],beta2v[j],beta3v[j],
                         kappa2v,gamma3v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    cat(l1,"-> ")
    
    # Stage 1D optimise over the kappa2(t)
    for(i in 1:m)
    {		 		 
      # cycle through the range of years
      dv=dtx[i,]	# actual deaths
      ev=etx[i,]	# exposure
      kappa2v[i]=llmaxM2D(beta1v,beta2v,beta3v,
                          kappa2v[i],gamma3v[(n+i-1):i],dv,ev,wv=wa[i,])
    }
    
    # Now apply the constraints
    fac21=mean(kappa2v)
    fac22=sum(beta2v)
    kappa2v=fac22*(kappa2v-fac21)    # ensures that the kappas sum to 0
    beta2v=beta2v/fac22		     # ensures that the beta2's sum to 1
    beta1v=beta1v+beta2v*fac22*fac21 # => beta1 needs to be adjusted to compensate
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    cat(l1," ->")
    
    # Stage 1A optimise over the beta1(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      wv=1	    # can be set to a vector of weights
      # to e.g. exclude duff years
      wv=wa[,j]
      s1=sum(wv*dtx[,j])
      s2=sum(wv*etx[,j]*exp(beta2v[j]*kappa2v+beta3v[j]*gamma3v[(n+1-j):(n+m-j)]))
      beta1v[j]=log(s1)-log(s2)
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    cat(l1,"\n")
    
  }		 # end while loop
  
  # calculate number of parameters and deduct 4 for the number of constraints
  npar=length(beta1v)+length(beta2v)+length(kappa2v)-2
  
  # Calculate the BIC
  BIC=l1-0.5*log(sum(wa))*npar
  
  list(beta1=beta1v,beta2=beta2v,beta3=beta3v,
       kappa2=kappa2v,gamma3=gamma3v,x=xv,y=yv,cy=cy,
       wa=wa,epsilon=epsilon,mhat=mhat,ll=l1,BIC=BIC,npar=npar,mtxLastYear = mtx[m,])		 
}

Df <- read.table('Df.csv')

ages <- 0:110  # Define ages up to the maximum supported age
years <- 1970:2019

lx_subset <- Df$lx[Df$Year %in% years & Df$Age %in% ages]
# Assuming you have already defined years and ages vectors
# Create a matrix from lx_subset with dimensions m*n
lx_matrix <- matrix(lx_subset, nrow = length(years), ncol = length(ages), byrow = TRUE)
# Optionally, you can set row and column names
rownames(lx_matrix) <- years
colnames(lx_matrix) <- ages
# Now lx_matrix is a matrix with dimensions m*n

dx_subset <- Df$dx[Df$Year %in% years & Df$Age %in% ages]
# Assuming you have already defined years and ages vectors
# Create a matrix from dx_subset with dimensions m*n
dx_matrix <- matrix(dx_subset, nrow = length(years), ncol = length(ages), byrow = TRUE)
# Optionally, you can set row and column names
rownames(dx_matrix) <- years
colnames(dx_matrix) <- ages

LCfit701 <- fit701(ages, years, lx_matrix, dx_matrix, matrix(1, length(years), length(ages)))

## ------------------------------------------------------------------------------------
names(LCfit701)

## ------------------------------------------------------------------------------------
library(ggplot2)

# Fit a Poisson regression model
model <- glm(dx ~ Year + Age, data = Df, family = "poisson")

# Predictions
Df$predicted_dx <- predict(model, type = "response")

# Plotting the results
p <- ggplot(Df, aes(x = Age, y = dx)) +
  geom_point() +
  geom_line(aes(y = predicted_dx), color = "red") +
  labs(title = "Poisson Regression Model", x = "Age", y = "dx")

print(p)
