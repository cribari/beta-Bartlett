##########################################################################################
#
#  PROGRAM: application1.r								                  
#  
#  USAGE: Computation of the corrected  test statistics for the data on 
#         the fifty countries with the largest prevalence of atheists 
#         from  Cribari-Neto and Souza (2013). 
#
#  NULL HYPOTHESIS: H_0: beta_4 = 0                                   								                              
# 
#  MODEL: g_1(mu_i)   = beta_1 + beta_2 x_{i2} + beta_3 x_{i3} + beta_4 x_{i4}  
#                     + beta_5 x_{i5} + beta_6 x_{i6} + beta_7 x_{i7} + beta_8 x_{i8}  
#         g_2(\phi_i) = delta_1 + delta_2 x_{i2}                                    	  
#
#  AUTHOR: Cristina Guedes 									  																		  
#########################################################################################

# R packages  required
require(betareg)

# Function to Bartlett's correction factor 
source("beta_bartlett.r")

# Reading data 
data <- read.table("atheists-data.txt", header=TRUE)
summary(data)

# Variables used in the analysis 
y      <- data$god/100
IQ     <- data$IQ
ppp1   <- data$ppp/1000
IQ2    <- data$IQ^2
tradeL <- log(data$trade)
int1   <- ppp1*tradeL
life   <- data$life
int2   <- life*tradeL

# Sample size
n      <- length(y)

# Link function: logit (mean)
glik   <- function(x){return(log(x/(1-x))) }
glin   <- function(x){return(1/(x*(1-x)))}
g2lin  <- function(x){return( (2*x-1)/((x^2)*(1-x)^2) )}
g3lin  <- function(x){return( (2 - 8*x + 12*x^2 - 6*x^3)/((x^3)*(1-x)^4) )}

# Link function: square root (precision)
hlik   <- function(x){return(sqrt(x))}
hlin   <- function(x){return( 1/(2*x^(1/2)) )}
h2lin  <- function(x){return(-( 1/(4*x^(3/2)) ))}
h3lin  <- function(x){return(( 3/(8*x^(5/2)) ))}

# Matrices of covariates (irrestricted) 
mX    <- cbind(1,IQ,IQ2,life,tradeL,ppp1,int1,int2) 
mZ    <- cbind(1,IQ)

# Matrices of covariates (restricted) 
mXr   <- cbind(1,IQ,IQ2,tradeL,ppp1,int1,int2)
mZ    <- cbind(1,IQ) 

# Covariate indices under H_0  
indice_X_H0 <- c(1,2,3,5,6,7,8) 
indice_Z_H0 <- c(9,10) 

# Number of parameters 
p   <- ncol(mX)
q   <- ncol(mZ)
k   <- p+q 

# Number of parameters of interest 
r <- ncol(mX)+ncol(mZ)-length(indice_X_H0)-length(indice_Z_H0) 

# Transformed variable
y_star  <- qlogis(y)
y_dag   <- log(1-y)

#################################################################################
##########                 BETA REGRESSION ESTIMATION                ############  
################################################################################# 

########## Unrestricted

# Initial values 
cf      <- coef(lm(y_star~mX[,2]+mX[,3]+mX[,4]+mX[,5]+mX[,6]+mX[,7]+mX[,8]))
etaols  <- mX%*%cf
muols   <- exp(etaols)/(1+exp(etaols))
varols  <- as.numeric(crossprod(y_star-etaols)) / ((n-p) * ((1/(muols*(1-muols)))^2))
phiols  <- ( muols * (1.0-muols) / varols )
cf1     <-  coef(lm(sqrt(phiols)~mZ[,2]))

# Model estimation
estima_full <- betareg(y~mX[,2]+mX[,3]+mX[,4]+mX[,5]+mX[,6]+mX[,7]+mX[,8]  | mZ[,2], 
                       start=c(cf,cf1), link="logit", link.phi="sqrt", data = data) 

########## Restricted

# Initial values 
cfr      <- coef(lm(y_star~mX[,2]+mX[,3]+mX[,5]+mX[,6]+mX[,7]+mX[,8]))
etaolsr  <- mXr%*%cfr
muolsr   <- exp(etaolsr)/(1+exp(etaolsr))
varolsr  <- as.numeric(crossprod(y_star-etaolsr)) / ((n-p-r) * ((1/(muolsr*(1-muolsr)))^2))
phiolsr  <- ( muolsr * (1.0-muolsr) / varolsr )
cf1r     <-  coef(lm(sqrt(phiolsr)~mZ[,2]))

# Model estimation 
estima_null <- betareg(y~mX[,2]+mX[,3]+mX[,5]+mX[,6]+mX[,7]+mX[,8] | mZ[,2], start=c(cfr,cf1r), 
                       link="logit", link.phi="sqrt", data = data) 

# Quantities used for calculation 
mu_hat_full     <- estima_full$fitted.values 
zeta_hat_full   <- estima_full$coefficients$precision %*% t(mZ)
phi_hat_full    <- (zeta_hat_full^2)[1:n] 
mcov_full       <- estima_full$vcov 
loglik_full     <- estima_full$loglik 

mu_hat_null     <- estima_null$fitted.values 
zeta_hat_null   <- estima_null$coefficients$precision %*% t(mZ) 
phi_hat_null    <- (zeta_hat_null^2)[1:n] 
mcov_null       <- estima_null$vcov 
loglik_null     <- estima_null$loglik 

# Likelihood ratio test statistics
LR <- 2*(loglik_full-loglik_null) # omega

# Model selection criteria 
AIC  <- AIC(estima_full)
BIC  <- BIC(estima_full)
AICc <- -2.0*loglik_full + 2.0*(k)*(n/(n-(k)-1))

# Pseudo R-squared
pseudoR2 <- estima_full$pseudo.r.squared

#################################################################################
##########        SKOVGAARD ADJUSTMENT (BETA REGRESSION)             ############  
#################################################################################
 
Y_star <- diag(y_star)
Y_dag  <- diag(y_dag)
iota   <- rep(1,n)
In     <- diag(rep(1,n)) 

########## Quantities under H_1
mu_star_full <- psigamma(mu_hat_full*phi_hat_full,deriv=0) - psigamma((1-mu_hat_full)*phi_hat_full,deriv=0)
mu_dag_full  <- psigamma((1-mu_hat_full)*phi_hat_full,deriv=0) - psigamma(phi_hat_full,deriv=0)
v_star_full  <- psigamma(mu_hat_full*phi_hat_full,deriv=1) + psigamma((1-mu_hat_full)*phi_hat_full,deriv=1) 
v_dag_full   <- psigamma((1-mu_hat_full)*phi_hat_full,deriv=1) - psigamma(phi_hat_full,deriv=1)
c_full       <- -psigamma((1-mu_hat_full)*phi_hat_full,deriv=1)  
M_full       <- diag(mu_hat_full)
M_star_full  <- diag(mu_star_full)
M_dag_full   <- diag(mu_dag_full)
Phi_full     <- diag(phi_hat_full)
T_full       <- diag(mu_hat_full * (1-mu_hat_full))
H_full       <- diag((2*phi_hat_full^(1/2)))

# Score function
U_b_full    <- t(mX) %*% Phi_full %*% T_full %*% (y_star - mu_star_full)
U_g_full    <- t(mZ) %*% H_full %*% (M_full %*% (y_star - mu_star_full) + (y_dag - mu_dag_full))
U_hat       <- rbind(U_b_full,U_g_full) # vetor escore

C_full      <- diag(c_full)
V_star_full <- diag(v_star_full)
V_dag_full  <- diag(v_dag_full)
S_full      <- diag((2*mu_hat_full - 1)/ ((mu_hat_full*(1-mu_hat_full))^2))
Q_full      <- diag(-( 1/(4*phi_hat_full^(3/2)) ))

# Fisher's information
I_hat      <- rbind(cbind( t(mX)%*%Phi_full%*%T_full%*%V_star_full%*%T_full%*%Phi_full%*%mX,
                     t(mX)%*%Phi_full%*%T_full%*%(M_full%*%V_star_full+C_full)%*%H_full%*%mZ),
              cbind( t(t(mX)%*%Phi_full%*%T_full%*%(M_full%*%V_star_full+C_full)%*%H_full%*%mZ),
                     t(mZ)%*%H_full%*%((M_full%*%V_star_full%*%M_full)+(2*M_full)%*%C_full+V_dag_full)%*%H_full%*%mZ))                 
I_hat_1    <- solve(I_hat) 

# Observed information
J_hat      <- rbind(cbind( t(mX)%*%(Phi_full%*%T_full%*%V_star_full+
                                S_full%*%(T_full%*%T_full)%*%(Y_star-M_star_full))%*%T_full%*%Phi_full%*%mX,
                     -t(mX)%*%((Y_star-M_star_full)-Phi_full%*%(M_full%*%V_star_full+C_full))%*%T_full%*%H_full%*%mZ),
              cbind( t(-t(mX)%*%((Y_star-M_star_full)-Phi_full%*%(M_full%*%V_star_full+C_full))%*%T_full%*%H_full%*%mZ),
                     t(mZ)%*%(H_full%*%((M_full%*%V_star_full%*%M_full)+(2*M_full)%*%C_full+V_dag_full)+
                    (M_full%*%(Y_star-M_star_full)+(Y_dag-M_dag_full))%*%Q_full%*%(H_full%*%H_full))%*%H_full%*%mZ))
J_hat_1    <- solve(J_hat)  


########## Quantities under H_0
mu_star_null <- psigamma(mu_hat_null*phi_hat_null,deriv=0) - psigamma((1-mu_hat_null)*phi_hat_null,deriv=0)
mu_dag_null  <- psigamma((1-mu_hat_null)*phi_hat_null,deriv=0) - psigamma(phi_hat_null,deriv=0)
v_star_null  <- psigamma(mu_hat_null*phi_hat_null,deriv=1) + psigamma((1-mu_hat_null)*phi_hat_null,deriv=1) 
v_dag_null   <- psigamma((1-mu_hat_null)*phi_hat_null,deriv=1) - psigamma(phi_hat_null,deriv=1)
c_null       <- -psigamma((1-mu_hat_null)*phi_hat_null,deriv=1) 
M_null       <- diag(mu_hat_null)
M_star_null  <- diag(mu_star_null)
M_dag_null   <- diag(mu_dag_null)
Phi_null     <- diag(phi_hat_null)
T_null       <- diag(mu_hat_null*(1-mu_hat_null))
H_null       <- diag((2*phi_hat_null^(1/2)))

# Score function
U_b_null    <- t(mX) %*% Phi_null %*% T_null %*% (y_star - mu_star_null)
U_g_null    <- t(mZ) %*% H_null %*% (M_null %*% (y_star - mu_star_null) + (y_dag - mu_dag_null))
U_til       <- rbind(U_b_null,U_g_null)

C_null      <- diag(c_null)
V_star_null <- diag(v_star_null)
V_dag_null  <- diag(v_dag_null)
S_null      <- diag((2*mu_hat_null - 1) / (mu_hat_null*(1-mu_hat_null))^2)
Q_null      <- diag(-( 1/(4*phi_hat_null^(3/2)) ))

# Fisher's information
I_til       <- rbind(cbind( t(mX)%*%Phi_null%*%T_null%*%V_star_null%*%T_null%*%Phi_null%*%mX,
                     t(mX)%*%Phi_null%*%T_null%*%(M_null%*%V_star_null+C_null)%*%H_null%*%mZ),
              cbind( t(t(mX)%*%Phi_null%*%T_null%*%(M_null%*%V_star_null+C_null)%*%H_null%*%mZ),
                     t(mZ)%*%H_null%*%(((M_null^2)%*%V_star_null)+(2*M_null)%*%C_null+V_dag_null)%*%H_null%*%mZ))
I_til_1     <- solve(I_til)

# Observed information
J_til       <- rbind(cbind( t(mX)%*%(Phi_null%*%T_null%*%V_star_null+
                                S_null%*%(T_null^2)%*%(Y_star-M_star_null))%*%T_null%*%Phi_null%*%mX,
                     -t(mX)%*%((Y_star-M_star_null)-Phi_null%*%(M_null%*%V_star_null+C_null))%*%T_null%*%H_null%*%mZ),
              cbind( t(-t(mX)%*%((Y_star-M_star_null)-Phi_null%*%(M_null%*%V_star_null+C_null))%*%T_null%*%H_null%*%mZ),
                     t(mZ)%*%(H_null%*%(((M_null^2)%*%V_star_null)+(2*M_null)%*%C_null+V_dag_null)+
                      (M_null%*%(Y_star-M_star_null)+(Y_dag-M_dag_null))%*%Q_null%*%(H_null^2))%*%H_null%*%mZ))
J_til_1     <- solve(J_til)  

# Upsilon_bar
Upsi_bar    <- rbind(cbind( t(mX)%*%Phi_full%*%T_full%*%V_star_full%*%T_null%*%Phi_null%*%mX ,
                        t(mX)%*%Phi_full%*%T_full%*%(V_star_full%*%M_null+C_full)%*%H_null%*%mZ ),
                 cbind( t(mZ)%*%H_full%*%(M_full%*%V_star_full+C_full)%*%T_null%*%Phi_null%*%mX,
                        t(mZ)%*%H_full%*%(M_full%*%V_star_full%*%M_null+(M_full+M_null)%*%C_full+V_dag_full)%*%H_null%*%mZ ))
Upsi_bar_1  <- solve(Upsi_bar)  # iversa de upsilon

# q_bar
q_bar       <- rbind( t(mX)%*%Phi_full%*%T_full%*%(V_star_full%*%(Phi_full%*%M_full-Phi_null%*%M_null)+C_full%*%(Phi_full-Phi_null))%*%iota,                  
               t(mZ)%*%H_full%*%((M_full%*%V_star_full+C_full)%*%(Phi_full%*%M_full-Phi_null%*%M_null)+(M_full%*%C_full+V_dag_full)%*%(Phi_full-Phi_null))%*%iota)

# xi 
csi        <- (((abs(det(I_til))*abs(det(I_hat))*abs(det(J_til[c(indice_X_H0,indice_Z_H0),c(indice_X_H0,indice_Z_H0)])))^(1/2))/
              (abs(det(Upsi_bar))*abs(det((I_til%*%Upsi_bar_1%*%J_hat%*%I_hat_1%*%Upsi_bar)[c(indice_X_H0,indice_Z_H0),c(indice_X_H0,indice_Z_H0)]))^(1/2)))*
              ((abs((t(U_til)%*%Upsi_bar_1%*%I_hat%*%J_hat_1%*%Upsi_bar%*%I_til_1%*%U_til))^(r/2))/
               abs(((LR^((r/2)-1))*t(U_til)%*%Upsi_bar_1%*%q_bar ))) 

# Skovgaard's modified likelihood ratio test statistic
LR_skov1   <- LR-2*log(csi[1])             # omega_{a1}
LR_skov2   <- LR*(1-(1/LR)*log(csi[1]))^2  # omega_{a2}

#################################################################################
######  BARTLETT CORRECTION TO THE LIKELIHOOD RATIO TEST (BETA REGRESSION) ######  
#################################################################################

epsilon_k   <- beta_bartlett(mu_hat_full,phi_hat_full ,mcov_full,mX,mZ)
epsilon_k_q <- beta_bartlett(mu_hat_null,phi_hat_null,mcov_null,mXr,mZ)  

B    <- epsilon_k-epsilon_k_q

# Bartlett's correction factor
c    <- 1+B/r   

LRs1 <- LR/c         # omega_{b1}
LRs2 <- LR*exp(-B/r) # omega_{b2}
LRs3 <- LR*(1-B/r)   # omega_{b3}

#################################################################################
#########                           RESULTS                          ############  
#################################################################################

# Program details 
program  <- rbind(n,"atheists-data.txt", estima_full$converged, estima_null$converged)
rownames(program) <- c("Sample size:", "Data:","Convergence (unrestricted):", 
                "Convergence (restricted):" )
colnames(program) <- c("")

# Parameter estimates and asymptotic standard errors:
Estim    <- rbind(cbind(estima_full$coefficients$mean),cbind(estima_full$coefficients$precision))
EP       <- cbind(diag((estima_full$vcov)^(1/2)) )
pv_ztest      <- 2.0*(1.0-pnorm(abs(Estim/EP))) 
Est      <- cbind(round(Estim,4),round(EP,4),pv_ztest)
colnames(Est) <- c("Estimate", "Stand. errors", "p-value(z-test)")
rownames(Est) <- c("beta_1", "beta_2","beta_3","beta_4","beta_5","beta_6","beta_7",
                   "beta_8","delta_1","delta_2")

# Value of test statistics and p-value (omega, omega_{b3} and omega_{a1})
resul   <- rbind(c(LR  ,pchisq(LR,r, lower.tail = FALSE)),
               c(LRs3,pchisq(LRs3,r, lower.tail = FALSE)),
               c(LR_skov1,pchisq(LR_skov1,r, lower.tail = FALSE))
              )
static_pv <- rbind(round(resul,4))
rownames(static_pv)<-c("w","wb3","wa1")
colnames(static_pv)<-c("Value","p-value")

# Fitted model quality measures
med      <- rbind(round(c(AIC, AICc, BIC,pseudoR2),4))
colnames(med)<-c("AIC","AICc","BIC","PseudoR2")
rownames(med)<-c("")

# Summary of results 
results  <- list(Program=program, Parameter_estimates=round(Est,4) ,Test_statistic = static_pv, 
                 Quality_measures= med)
results


                                                



