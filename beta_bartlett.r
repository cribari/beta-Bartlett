############################################################################################
#
#  PROGRAM:   beta_bartlett.r
#
#  USAGE:    Given some input arguments for the function, which come from the beta 
#            regression model with adjusted variable dispersion, the function returns the 
#            epsilon_k of Bartlett's correction. This function uses the general matrix 
#            expression to obtain the bartlett correction factor presented in Cordeiro (1993). 
#
#  NOTE:     In a simulation study, this function must be called within another file with 
#            the implementation of interest. 
#
#  AUTORA:    Cristina Guedes 
#
#############################################################################################

# Auxiliary function (zero counter)
manyzeros<-function(r)
{
  zeros <- 0
  for( i in 1:length(r))
  {
    if(r[i]==0) zeros<-zeros+1
  } 
  return(zeros)
}

# Function that determines epsilon_k of Bartlett's correction factor 
beta_bartlett <- function(mut,phit,K_1,X,Z)
{  
  ################################################
  # Function input arguments
  # mut - vector mu 
  # phit - vector phi
  # K_1 - inverse of the information matrix
  # X - matrix (regressors of mean)
  # Z - matrix (regressors of precision)
  ################################################ 
  
  # Link function: logit (mean)
  glik  <- function(x){return(log(x/(1-x))) }
  glin  <- function(x){return(1/(x*(1-x)))}
  g2lin <- function(x){return( (2*x-1)/((x^2)*(1-x)^2) )}
  g3lin <- function(x){return( (2 - 8*x + 12*x^2 - 6*x^3)/((x^3)*(1-x)^4))}                                                              
  
  # Link function: square root (precision)
  hlik  <- function(x){return(sqrt(x))}
  hlin  <- function(x){return( 1/(2*x^(1/2)) )}
  h2lin <- function(x){return(-( 1/(4*x^(3/2)) ))}
  h3lin <- function(x){return(( 3/(8*x^(5/2)) ))}
  
  # Number of parameters 
  pf       <- ncol(X) # number of parameters of the average submodel 
  qf       <- ncol(Z) # number of parameters of the precision submodel 
  kf       <- pf+qf 
  indice_k <- seq(1,kf) 
  
  # General matrix of regressors 
  XZ=cbind(X,Z) 
  
  # Auxiliary function (deltas counter)
  contador <- function(r)
  {
    zeros <- 0
    for( i in 1:length(r))
    {
      if(as.numeric(r[i]>pf)==1) zeros<-zeros+1
    } 
    return(zeros)
  }
  
  # Required quantities 
  wt           <- psigamma(mut*phit,deriv=1) + psigamma((1-mut)*phit,deriv=1)
  mt           <- psigamma(mut*phit,deriv=2) - psigamma((1-mut)*phit,deriv=2)
  DmuDmueta    <- -g2lin(mut)/(glin(mut)^2)
  DmuDeta      <- 1/glin(mut)
  DphiDphizeta <-  -h2lin(phit)/(hlin(phit)^2) 
  DphiDzeta    <- 1/hlin(phit) 
  at           <- 3*DmuDmueta*(DmuDeta^2)
  tt           <- 3*DphiDphizeta*(DphiDzeta^2)
  
  dt           <- ((1-mut)^2)*psigamma((1-mut)*phit,deriv=1)+ (mut^2)*psigamma(mut*phit,deriv=1)-psigamma(phit,deriv=1)
  st           <- ((1-mut)^3)*psigamma((1-mut)*phit,deriv=2)+ (mut^3)*psigamma(mut*phit,deriv=2)-psigamma(phit,deriv=2)
  ct           <- phit*(mut*wt - psigamma((1-mut)*phit,deriv=1))
  DwDphi       <- mut*psigamma(mut*phit,deriv=2) + (1-mut)*psigamma((1-mut)*phit,deriv=2)
  DmustarDmu   <- phit*wt
  DmudagDmu    <- -phit*psigamma((1-mut)*phit,deriv=1)
  DmustarDphi  <- ct/phit
  DmudagDphi   <- (1-mut)*psigamma((1-mut)*phit,deriv=1)-psigamma(phit,deriv=1)
  D2mustarDphi <- (mut^2)*psigamma(mut*phit,deriv=2)-((1-mut)^2)*psigamma((1-mut)*phit,deriv=2)
  D3mustarDphi <- (mut^3)*psigamma(mut*phit,deriv=3)-((1-mut)^3)*psigamma((1-mut)*phit,deriv=3)
  ut           <- -phit*(2*wt+phit*(DwDphi))
  rt           <- (2*(DmustarDphi) + phit*D2mustarDphi)*DmuDeta
  zt           <- DmustarDphi + phit*D2mustarDphi
  Dmueta3      <- -3*g2lin(mut)/(glin(mut)^4)
  Dmueta2      <- -2*g2lin(mut)/(glin(mut)^3)
  D2mueta      <- (-g3lin(mut)*glin(mut)+2*(g2lin(mut)^2))/(glin(mut)^3)
  Dphizeta3    <-  -3*h2lin(phit)/(hlin(phit)^4) #
  Dphizeta2    <-  -2*h2lin(phit)/(hlin(phit)^3) #
  D2phizeta    <-  (-h3lin(phit)*hlin(phit)+2*(h2lin(phit)^2))/(hlin(phit)^3) 
  bt           <- DmuDeta*(D2mueta*DmuDeta + (DmuDmueta^2))
  vt           <- DphiDzeta*(D2phizeta*DphiDzeta + (DphiDphizeta^2))
  DaDmu        <- 3*DmuDeta*(D2mueta*DmuDeta+2*(DmuDmueta^2))
  DtDphi       <- 3*DphiDzeta*(D2phizeta*DphiDzeta+2*(DphiDphizeta^2))
  DmDphi       <- mut*psigamma(mut*phit,deriv=3) - (1-mut)*psigamma((1-mut)*phit,deriv=3)
  D2wDphi      <- (mut^2)*psigamma(mut*phit,deriv=3) + ((1-mut)^2)*psigamma((1-mut)*phit,deriv=3)
  DrDmu        <- (2*DmustarDphi + phit*D2mustarDphi)*DmuDmueta + (2*wt+4*phit*DwDphi+(phit^2)*D2wDphi)*DmuDeta
  DsDmu        <- 3*D2mustarDphi + phit*D3mustarDphi
  DsDphi       <- (mut^4)*psigamma(mut*phit,deriv=3) + ((1-mut)^4)*psigamma((1-mut)*phit,deriv=3) - psigamma(phit,deriv=3)
  DuDphi       <- -2*wt - 4*phit*DwDphi-(phit^2)*D2wDphi
  DrDphi       <- DsDmu*DmuDeta
  DcDmu        <- phit*(wt + phit*DwDphi)
  DuDmu        <- -(phit^2)*(3*mt + phit*DmDphi)
  DzDphi       <- 2*D2mustarDphi + phit*D3mustarDphi
  DzDmu        <- wt + 3*phit*DwDphi + (phit^2)*D2wDphi
  DmDmu        <- phit*psigamma(mut*phit,deriv=3) + phit*psigamma((1-mut)*phit,deriv=3)
  DmuDwphi     <- mt + phit*DmDphi
  DwDmu        <- phit*mt
  
  #########  Cumulants of the varying precision beta regression model 
  
  k4 <- function(r,s,t,u) # k_rstu
  {
    vetor_theta   <- c(r,s,t,u) # theta indices vector 
    ordem_indices <- order(vetor_theta,decreasing=FALSE)
    
    ndelta <- contador(vetor_theta) 
    
    if(ndelta == 0) # k_rstu
    {
      kappa <- sum(-phit*(phit^2*(mt*Dmueta3 + DmDmu*(DmuDeta)^3 ) + 
                           phit*(wt*(DaDmu + bt) + DwDmu*at))*DmuDeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
    }
    if(ndelta == 1) # k_rstR
    {
      kappa <- sum(-phit*(phit*(3*mt+phit*DmDphi)*DmuDeta^3 + at*(2*wt+phit*DwDphi) + 
                           bt*ct/phit)*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
    }
    if(ndelta == 2) # k_rsRS
    {
      kappa <- sum((-ct*DmuDmueta*DphiDphizeta*DphiDzeta - DrDmu*(DphiDzeta^2) - DmustarDphi*DmuDeta*(DphiDzeta^2) + ut*DmuDeta*DphiDzeta*DphiDphizeta)*DmuDeta*
                    XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
    }
    if(ndelta == 3) # k_rRST
    {
      kappa<- sum((-ct*DmuDeta*vt -rt*Dphizeta3 - DsDmu*(DphiDzeta^3)*DmuDeta)*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
      
    }
    if(ndelta == 4) # k_RSTU
    {
      kappa <- sum((-DsDphi*(DphiDzeta^3) - 2*st*Dphizeta3 - dt*(DtDphi + vt))*
                    DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
    }
    return(kappa)
  }
  
  k3 <- function(r,s,t) # k_rst
  {
    vetor_theta   <- c(r,s,t) 
    ordem_indices <- order(vetor_theta,decreasing=FALSE) 
    
    ndelta <- contador(vetor_theta) # quantos deltas's 
    
    if(ndelta == 0) # k_rst
    {
      kappa <- sum(-phit^2*(phit*mt*(DmuDeta^3) + wt*at)*
                    XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
    }
    if(ndelta == 1) # k_rsR
    {
      kappa <- sum((ut*DmuDeta - ct*DmuDmueta)*DmuDeta*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
    }
    if(ndelta == 2) # k_rRS
    {
      kappa <- sum((-rt*DphiDzeta - ct*DphiDphizeta*DmuDeta)*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
    }
    if(ndelta == 3) # k_RST
    {
      kappa <- sum((-st*(DphiDzeta^3) - dt*Dphizeta2*DphiDzeta - dt*DphiDphizeta*(DphiDzeta^2))*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
      
    }
    return(kappa)
  }
  
  k31 <- function(r,s,t,U) # k_rst(u)
  {
    vetor_theta   <- c(r,s,t,U) 
    ordem_indices <- order(vetor_theta,decreasing=FALSE) 
    
    ndelta <- contador(vetor_theta) 
    
    if(as.numeric(U>pf)==1) # k_rst(R)
    {
      if(ndelta == 1) # k_rst(R)
      {
        kappa <- sum(-phit*((DmuDeta^3)*(3*phit*mt + phit^2 *DmDphi) + at*(2*wt + phit*DwDphi))*DphiDzeta*
                      XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
      }
      if(ndelta == 2) # k_rsR(S)
      {
        kappa <- sum(((DphiDzeta*DuDphi + ut*DphiDphizeta)*(DmuDeta^2) - (DphiDzeta*zt + ct*DphiDphizeta)*DmuDmueta*DmuDeta)*DphiDzeta*
                      XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
      }
      if(ndelta == 3) # k_rRS(T)
      { 
        
        kappa <- sum((-DrDphi*(DphiDzeta^3) - rt*Dphizeta2*DphiDzeta - zt*DmuDeta*DphiDphizeta*(DphiDzeta^2) - ct*DmuDeta*vt)* XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
      }
      if(ndelta == 4) # k_ RST(U)
      {
        kappa <- sum((-DsDphi*(DphiDzeta^3) - 2*st*Dphizeta3 - dt*DtDphi)*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
      }    
    }else{ # k_rst(u), if U>0
      
      if(ndelta == 0) # k_rst(u)
      {
        kappa <- sum(-phit^2 *(phit*(mt*(Dmueta3 + at) + (DmuDeta^3)*DmDmu) + wt*DaDmu)*DmuDeta*
                      XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
      }
      if(ndelta == 1) # k_rsR(t)
      {  
        kappa <- sum((DuDmu*(DmuDeta^3) + 2*ut*DmuDmueta*(DmuDeta^2) - DcDmu*DmuDmueta*(DmuDeta^2)-ct*bt)*DphiDzeta*
                      XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])

      }
      if(ndelta == 2) # k_rRS(s)
      {
        kappa <- sum((-DrDmu*DphiDzeta^2 - (DcDmu*DmuDeta + ct*DmuDmueta)*DphiDphizeta*DphiDzeta)*DmuDeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
      }
      if(ndelta == 3) # k_RST(r)
      {
        kappa <- sum((-DsDmu*(DphiDzeta^3) - (2*DmustarDphi + phit*D2mustarDphi)*Dphizeta3)*DmuDeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]
                    *XZ[,vetor_theta[ordem_indices[4]]])
        
      } 
    }
    return(kappa)
  }
  
  k22<-function(r,s,T,U) # k_rs(tu)
  {
    vetor_theta   <- c(r,s,T,U) 
    ordem_indices <- order(vetor_theta,decreasing=FALSE) 
    
    ndelta<-contador(vetor_theta)  
    
    if(as.numeric(U>pf)==1)
    {
      if(as.numeric(T>pf)==1)# k_rs(RS)
      {
        if(ndelta == 2) # k_rs(RS)
        {
          kappa <- sum(((DuDphi*DphiDzeta + ut*DphiDphizeta)*(DmuDeta^2))*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]]) 
        }
        if(ndelta == 3) # k_rR(ST)
        {
          kappa <- sum((-DzDphi*(DphiDzeta^3)*DmuDeta - zt*Dphizeta3*DmuDeta  - ct*vt*DmuDeta)* 
                        XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
  
        }
        if(ndelta == 4) # k_RS(TU)
        {
          kappa <- sum((-DsDphi*(DphiDzeta^3) - st*Dphizeta3 - (2/3)*st*tt - (2/3)*dt*DtDphi)*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]
                      *XZ[,vetor_theta[ordem_indices[4]]])
        }
      }else{ # k_rs(tR)
        if(ndelta == 1) # k_rs(tR) 
        { 
          kappa <- sum((DuDmu*DmuDeta^3 + (2/3)*ut*at)*DphiDzeta*
                        XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
        }
        if(ndelta == 2) # k_rR(sS)
        {
          kappa <- sum((-(wt + 3*phit*DwDphi+ (phit^2)*D2wDphi)*DmuDeta^2*DphiDzeta - DcDmu*DmuDeta^2*DphiDphizeta - zt*DmuDmueta*DphiDzeta*DmuDeta -
                         ct*DmuDmueta*DmuDeta*DphiDphizeta)*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
        }
        if(ndelta == 3) # k_RS(rT)
        {
          kappa <- -sum((DrDphi*DphiDzeta^2 + rt*Dphizeta2)*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]
                       *XZ[,vetor_theta[ordem_indices[4]]])
        }
      }
    }else{ # if U > 0
      if(as.numeric(T>pf)==1) # k_rs(Rt)
      {
        if(ndelta == 1) # k_rs(Rt)
        { 
          kappa <- sum((DuDmu*DmuDeta^2 + ut*Dmueta2)*(DmuDeta)*DphiDzeta*
                        XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])

        }
        if(ndelta == 2) # k_rR(Ss)
        {
          kappa <- -sum((DzDmu*DmuDeta*(DphiDzeta^2) + zt*DmuDmueta*(DphiDzeta^2) + DcDmu*DphiDphizeta*DphiDzeta*DmuDeta +
                          ct*DphiDphizeta*DphiDzeta*DmuDmueta)*DmuDeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]
                       *XZ[,vetor_theta[ordem_indices[4]]])
        }
        if(ndelta == 3) # k_RS(Tr)
        {
          kappa <- -sum((DrDphi*DphiDzeta^2 + rt*Dphizeta2)*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]
                       *XZ[,vetor_theta[ordem_indices[4]]])
        }
      }else{ # k_rs(tu)
        if(ndelta == 0) # k_rs(tu)
        {
          kappa <- sum(-(phit^2)*(phit*(mt*(Dmueta3 + (2/3)*at) + (DmuDeta^3)*DmDmu) + (2/3)*wt*DaDmu)*DmuDeta*
                        XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
        }
        if(ndelta == 1) # k_rR(st)
        {
          kappa <- sum(((-phit^2)*(mt + DmuDwphi)*(DmuDeta^3) - DcDmu*Dmueta3 - ct*bt)*DphiDzeta*
                        XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]*XZ[,vetor_theta[ordem_indices[4]]])
        }
        if(ndelta == 2) # k_RS(rs)
        {
          kappa <- -sum((DrDmu*DmuDeta*(DphiDzeta^2))*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]]
                       *XZ[,vetor_theta[ordem_indices[4]]])
        }
      }
    }
    return(kappa) 
  }
  
  k21 <- function(r,s,T) # k_rs(t)
  {
    vetor_theta<-c(r,s,T) 
    ordem_indices<-order(vetor_theta,decreasing=FALSE) 
    
    ndelta <- contador(vetor_theta) 
    
    if(as.numeric(T>pf)==1) # k_rs(R)
    {
      if(ndelta == 1) # k_rs(R)
      {
        kappa <- sum((ut*(DmuDeta^2)*DphiDzeta)*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
      }
      if(ndelta == 2) # k_rR(S)
      {
        kappa <- sum((-zt*DphiDzeta - ct*DphiDphizeta)*DmuDeta*DphiDzeta*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
      }
      if(ndelta == 3) # k_RS(T)
      {
        kappa <- -sum((st*(DphiDzeta^3) + (2/3)*dt*tt)*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
      }  
    }else{
      if(ndelta == 0) # k_rs(t)
      {
        kappa <- sum(-(phit^2)*(phit*mt*(DmuDeta^3) + (2/3)*wt*at)*
                      XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
      }
      if(ndelta == 1) # k_rR(s) 
      {
        kappa <- sum((-DcDmu*DmuDeta*DphiDzeta - ct*DmuDmueta*DphiDzeta)*DmuDeta*
                      XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
      }
      if(ndelta == 2) # k_RS(r)
      {
        kappa <- -sum(rt*(DphiDzeta^2)*XZ[,vetor_theta[ordem_indices[1]]]*XZ[,vetor_theta[ordem_indices[2]]]*XZ[,vetor_theta[ordem_indices[3]]])
      }
    }
    return(kappa) 
  }
  
  #########  The general adjustment factor in matrix form
  # matrix initialization
  mA <- array(rep(0,kf^4),dim=c(kf,kf,kf,kf))
  mP <- array(rep(0,kf^3),dim=c(kf,kf,kf))
  mQ <- mP
  
  # Matrix mA
  for(ti in indice_k)
  {
    
    for(ui in indice_k)
    {
      
      for(ri in indice_k)
      {
        
        for(si in indice_k)
        {
          
          mA[ti,ui,ri,si] <- k4(ri,si,ti,ui)/4 - k31(ri,si,ti,ui) + k22(ri,ti,si,ui)
        }
      }
    }
  }
  
  # Matrix mP
  for(ti in indice_k)
  {
    
    for(ri in indice_k)
    {
      
      for(si in indice_k)
      {
        
        mP[ti,ri,si]<- k3(ri,si,ti)
      }
    }
  }
  
  # Matrix mQ
  for(ui in indice_k)
  {
    
    for(ri in indice_k)
    {
      
      for(si in indice_k)
      {
        
        mQ[ui,ri,si]<- k21(si,ui,ri)	
      }
    }
  }
  
  # matrix initialization
  L<-matrix(rep(0,kf^2),ncol=kf)
  M1 <- L
  M2 <- L
  M3 <- L
  N1 <- L
  N2 <- L
  N3 <- L
  
  # Matrix L
  for(ri in indice_k)
  {
    for(si in indice_k)
    { 
      L[ri,si] <- sum(diag( K_1 %*% mA[ri,si,,] ))
    }
  }
  
  # Matrix Mi
  for(ri in indice_k)
  {
    for(si in indice_k)
    { 
      M1[ri,si] <- sum(diag( K_1%*%mP[ri,,]%*%K_1%*%mP[si,,] ))
    }
  }
  for(ri in indice_k)
  {
    for(si in indice_k)
    { 
      M2[ri,si] <- sum(diag( K_1%*%mP[ri,,]%*%K_1%*%t(mQ[si,,] )))
    }
  }
  for(ri in indice_k)
  {
    for(si in indice_k)
    { 
      M3[ri,si] <- sum(diag( K_1%*%mQ[ri,,]%*%K_1%*%mQ[si,,] ))
    }
  }
  
  # Matrix Ni
  for(ri in indice_k)
  {
    for(si in indice_k)
    { 
      N1[ri,si] <- (sum(diag( mP[ri,,]%*%K_1 )))*(sum(diag( mP[si,,]%*%K_1 )))
    }
  }
  for(ri in indice_k)
  {
    for(si in indice_k)
    { 
      N2[ri,si] <- (sum(diag( mP[ri,,]%*%K_1 )))*(sum(diag( mQ[si,,]%*%K_1 )))
    }
  }
  for(ri in indice_k)
  {
    for(si in indice_k)
    { 
      N3[ri,si] <- (sum(diag( mQ[ri,,]%*%K_1 )))*(sum(diag( mQ[si,,]%*%K_1 )))
    }
  }
  
  M <- -M1/6 + M2 - M3
  N <- -N1/4 + N2 - N3
  
  return(sum(diag( K_1%*%(L - M - N) )))  # epsilon_k of Bartlett's correction factor
  
} # end of function 

