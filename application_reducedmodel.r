##########################################################################################
#
#  PROGRAM: application_reducedmodel.r								                  
#
#  USAGE: Based on the testing inference, we arrive at the reduced model, which is our 
#         final model for the data on  the fifty countries with the largest prevalence of
#         atheists from  Cribari-Neto and Souza (2013).
#
#  REDUCED MODEL: g_1(mu_i)   = beta_1 + beta_2 x_{i2} + beta_3 x_{i3} + beta_4 x_{i6}  
#                     + beta_5 x_{i7} + beta_6 x_{i8}   
#                 g_2(\phi_i) = delta_1 + delta_2 x_{i2}                                    	  
#
#  AUTHOR: Cristina Guedes 									  																		  
#########################################################################################

# R packages  required
require(betareg)

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

# Matrices of covariates 
mX    <- cbind(1,IQ,IQ2,ppp1,int1,int2) 
mZ    <- cbind(1,IQ)

# Number of parameters 
p   <- ncol(mX)
q   <- ncol(mZ)
k   <- p+q

# Transformed variable
y_star  <- qlogis(y)

#################################################################################
##########                 BETA REGRESSION ESTIMATION                ############  
################################################################################# 

# Initial values 
cf      <- coef(lm(y_star~mX[,2]+mX[,3]+mX[,4]+mX[,5]+mX[,6]))
etaols  <- mX%*%cf
muols   <- exp(etaols)/(1+exp(etaols))
varols  <- as.numeric(crossprod(y_star-etaols)) / ((n-p) * ((1/(muols*(1-muols)))^2))
phiols  <- ( muols * (1.0-muols) / varols )
cf1     <-  coef(lm(sqrt(phiols)~mZ[,2]))

# Model estimation
estima_full <- betareg(y~mX[,2]+mX[,3]+mX[,4]+mX[,5]+mX[,6]  | mZ[,2], start=c(cf,cf1),
                       link="logit", link.phi="sqrt", data = data) 

# Model selection criteria 
AIC  <- AIC(estima_full)
BIC  <- BIC(estima_full)
AICc <- -2.0*(estima_full$loglik) + 2.0*(k)*(n/(n-(k)-1))

# Pseudo R-squared
pseudoR2 <- estima_full$pseudo.r.squared

#################################################################################
#########                           RESULTS                          ############  
#################################################################################

# Program details 
program  <- rbind(n,"atheists-data.txt", estima_full$converged)
rownames(program) <- c("Sample size:", "Data:","Convergence (unrestricted):" )
colnames(program) <- c("")

# Parameter estimates and asymptotic standard errors:
Estim    <- rbind(cbind(estima_full$coefficients$mean),cbind(estima_full$coefficients$precision))
EP       <- cbind(diag((estima_full$vcov)^(1/2)) )
pv_ztest      <- 2.0*(1.0-pnorm(abs(Estim/EP))) 
Est      <- cbind(round(Estim,4),round(EP,4),pv_ztest)
colnames(Est) <- c("Estimate", "Stand. errors", "p-value(z-test)")
rownames(Est) <- c("beta_1", "beta_2","beta_3","beta_4","beta_5","beta_6",
                   "delta_1","delta_2")

# Fitted model quality measures
med           <- rbind(round(c(AIC, AICc, BIC, pseudoR2),4))
colnames(med) <- c("AIC","AICc","BIC","PseudoR2")
rownames(med)<-c("")

# Summary of results 
results  <- list(Program=program, Parameter_estimates=round(Est,4) , Quality_measures= med)
results






