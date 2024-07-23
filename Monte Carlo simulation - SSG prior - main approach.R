# ==============================================================================
# ================================= Packages ===================================
# ==============================================================================
#Simulation SSG prior setup

library("wavethresh")
library("EbayesThresh")
library('coda')
library('bfw')
library('ggplot2')
library('latex2exp')
library('ggtext')

# ==============================================================================
# ====================== Creating the mixture behavior =========================
# ==============================================================================

#Defining the behavior for our dynamic mixture weight, \alpha_t.

sinusoidal <- function(q)
  cos(2*pi*(q+pi))/(2.5)+0.5
heavisinebe <- function(q){
  (6+4*sin(4*pi*q) - sign(q - 0.3) - sign(0.72 - q))/10}
stepwisebe <- function(q){
  ifelse(q<0.3, 0.2, ifelse(q<0.7, 0.8, 0.3))}
parabolicbe <- function(q){
  3*(q-0.5)^2 + 0.125}
staticbe <- function(q){
  rep(0.75, length(q))
}
bumpsbe <- function(q){
  pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65,  0.76, 0.78, 0.81)
  hgt <-  c( 4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  wth <- c(0.005,0.005,0.006,0.01, 0.01, 0.03, 0.01, 0.01,  0.005, 0.008,0.005)
  y5 <- rep(0, length(q))
  
  for(j in 1:length(pos)){
    y5 = y5 + hgt[j]/( 1 + abs((q - pos[j])/wth[j]))^4
  }
  y5 <- y5/5.1
  return(y5)
}
blocks <- function(q){
  pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65,  0.76, 0.78, 0.81)
  hgt <-  c( 4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
  y6 <- rep(0, length(q))
  for(j in 1:length(pos)){
    y6 = y6 + (1 + sign(q-pos[j]))*(hgt[j]/2)
  }
  y6 = (y6+2)/7.2
  return(y6)
}

#Since we are going to apply a Wavelet transform, we must choose a power of two 
#as our simulation size. Here we choose 1024.

serie <- seq_len(1024)/1024

# serie <- seq(1/1024, 1024/1024, 1/1024)

simulationsize <- length(serie)

#Creating our mixture weights, which will be stored in a vector called
#"alpha_behavior".

alpha_behavior <- sinusoidal(serie)
alpha_behavior <- heavisinebe(serie)
alpha_behavior <- stepwisebe(serie)
alpha_behavior <- parabolicbe(serie)
alpha_behavior <- staticbe(serie)
alpha_behavior <- bumpsbe(serie)
alpha_behavior <- blocks(serie)

# ============================ Global variables ================================
# Generating mixture parameters

#y_t = (1 - z_t)x_{1t}+z_t x_{2t}
#where x_{1t} \sim N(0,0.25); x_{2t} \sim N(2, 0.25)
#which means (\phi_1 = \phi_2 = 4)

truemean1 <- 0
truemean2 <- 2
trueprecision1 <- 4
trueprecision2 <- 4

#Gibbs sampler: Setting down the size of our chain, of our burn-in and our lag

burn <- 1000
lags <- 5
nchain <- 1000
BB <- burn + lags*nchain

#Family choice
family_choice <- "DaubExPhase" #"DaubLeAsymm", "DaubExPhase", "Coiflets"

#Number of vanishing moments
number_vm <- 10

# ==============================================================================
# ============================== Wavelet Function ==============================
# ==============================================================================

#Function to transform the wavelet coefficients back to the data domain
inverse_coeff <- function(coeff, vm, fc){
  zero <- rep(0, length(coeff))
  zero_wd <- wd(zero, filter.number = vm, family = fc)
  zero_wd$C[2^(log2(length(coeff))+1)-1] <- coeff[1]
  zero_wd$D <- coeff[2:length(coeff)]
  return(wr(zero_wd))    
}


# Function to calculate the posterior weight for non-null effect
wpost_normal <- function(w, xx, tau2){
  wpost <- (w*dnorm(xx, sd = sqrt(1+tau2)))/(w*dnorm(xx,
                                            sd = sqrt(1+tau2))+(1-w)*dnorm(xx))
  wpost[round(wpost, digits = 1)==0] <- 0
  wpost[1] <- 1
  return(wpost)
}

# ==============================================================================

zeta <- 1
rho <- 1
kappa <- 1
xi <- 100

alpha <- matrix(nrow=BB,ncol=simulationsize)
mu <- matrix(nrow=BB,ncol = 2)
phi <- matrix(nrow=BB,ncol = 2)

prior_alpha_gamma <- 0.01
prior_beta_gamma <- 0.01
sigma2_theta <- 100
sigma2_theta_post <- sigma2_theta/(sigma2_theta + 1)


N <- seq(burn+1, BB, lags)

num_rep <- 1000
estimated_alpha_median <- matrix(nrow=num_rep,ncol=simulationsize)
estimated_alpha_mean <- matrix(nrow=num_rep,ncol=simulationsize)
estimated_mu_median <- matrix(nrow=num_rep,ncol=2)
estimated_tau_median <- matrix(nrow=num_rep,ncol=2)

gc_after_a_number <- 1e2
# ==============================================================================
# ============================== Monte Carlo ===================================
# ==============================================================================

for(pp in 1:num_rep){
  
# ====================== Creating the mixture data =============================
  
  #Creating our latent variables Z_t, where Z_t|\alpha_t \sim Bern(\alpha_t)
  z <- rbinom(n = simulationsize, size = 1, prob = alpha_behavior)
  
  #y_t = (1 - z_t)x_{1t}+z_t x_{2t}
  
  x1 <- rnorm(n = simulationsize, mean = truemean1, sd = sqrt(1/trueprecision1))
  x2 <- rnorm(n = simulationsize, mean = truemean2, sd = sqrt(1/trueprecision2))
  
  y <- (1-z)*x1 + z*x2
  
  #Consider the following priors to the parameters \mu_k
  #\mu_1 \sim N(q_1, 10s^2)
  #\mu_2 \sim N(q_3, 10s^2)
  #\phi_k \sim \Gamma(0.01, 0.01) 
  prior_mean1 <- quantile(y, probs = 0.25)
  prior_mean2 <- quantile(y, probs = 0.75)
  prior_variance1 <- var(y)
  prior_variance2 <- var(y)
  
  #Defining our start values
  mu[1,]<-c(prior_mean1,prior_mean2)
  phi[1,]<- c(1/prior_variance1, 1/prior_variance2)
  coef_wavelet <- vector("double", simulationsize) 
  alpha[1,] <- rep(0.5, simulationsize)
  prob0start <- (1-alpha[1,])*dnorm(y, mean = mu[1,1], 
                                    sd = 1/sqrt(phi[1,1]))
  prob1start <- alpha[1,]*dnorm(y, mean = mu[1,2], 
                                sd = 1/sqrt(phi[1,2]))
  probvalstart <- prob1start/(prob0start + prob1start)
  
  pred_label <- rbinom(simulationsize, 1, prob = probvalstart)
  
  inv_coeff <- inverse_coeff(coef_wavelet, vm = number_vm, fc = family_choice) 
  lat_var <- rnorm(simulationsize, mean = inv_coeff, 1)
  pi_mix_wavelet <- rep_len(1, simulationsize)
  tau2_mix_gaussian <- rep_len(100, simulationsize)
  epsilon <- 1e-10  # Small constant to prevent division by zero
  mix_var <- rbinom(n = simulationsize, size = 1, prob = 0.1)
  ninho1 <- sum(mix_var)
  ninho0 <- simulationsize - ninho1
  
# ================================== Gibbs =====================================
  
  for(ii in 2:BB){
    #Creating some auxiliary variables for calculating the conditional posterior 
    n1 <- sum(pred_label)
    n0 <- simulationsize - n1
    y1 <- y[pred_label==1]
    y0 <- y[pred_label==0]
    sy1 <- sum(y1)
    sy0 <- sum(y0)
    
# ====================== Generating mean.lower =================================
    var_mu_lower <- 1/(n0*phi[ii-1,1]+ 1/prior_variance1)
    mean_mu_lower <- var_mu_lower *(sy0*phi[ii-1,1] +
                                      prior_mean1/prior_variance1)
    mu[ii,1] <- rnorm(n = 1, mean = mean_mu_lower, sd = sqrt(var_mu_lower))
    
# ====================== Generating prec.lower =================================  
    alpha_gamma_lower <- prior_alpha_gamma + n0/2
    beta_gamma_lower <- prior_beta_gamma + sum((y0 - mu[ii,1])^2)/2
    phi[ii,1] <- rgamma(1, shape = alpha_gamma_lower, rate = beta_gamma_lower)
    
# ====================== Generating mean.higher ================================  
    var_mu_higher <- 1/(n1*phi[ii-1,2]+ 1/prior_variance2)
    mean_mu_higher <- var_mu_higher *(sy1*phi[ii-1,2] +
                                        prior_mean2/prior_variance2)
    mu[ii,2] <- rnorm(n = 1, mean = mean_mu_higher, sd = sqrt(var_mu_higher))
    
# ====================== Generating prec.higher ================================ 
    alpha_gamma_higher <- prior_alpha_gamma + n1/2
    beta_gamma_higher <- prior_beta_gamma + sum((y1 - mu[ii,2])^2)/2
    phi[ii,2] <- rgamma(1, shape = alpha_gamma_higher, rate = beta_gamma_higher)
    
    
# ============================== Constraint ====================================
    #For dealing with label switching, we stablish that the pairs (\mu_k,\phi_k) 
    #are ordered under the constraint \mu_1 < \mu_2 and \phi_1 < \phi_2
    
    if(mu[ii,2] < mu[ii,1]){
      mu[ii,]<- mu[ii,2:1]
      phi[ii,]<- phi[ii,2:1]
    }
    
# ====================== Generating pred_label =================================
    prob0 <- (1-alpha[ii-1,])*dnorm(x = y, mean = mu[ii,1],
                                    sd = 1/sqrt(phi[ii,1]))
    prob1 <- alpha[ii-1,]*dnorm(x = y, mean = mu[ii,2],
                                sd = 1/sqrt(phi[ii,2]))
    probval <- prob1/(prob0 + prob1 + epsilon)
    
    pred_label <- rbinom(n = simulationsize, size = 1, prob = probval)
    
    
# ========================= Generating lat_var =================================
    U <- runif(n = simulationsize)
    upper_bound <- (pnorm(q = 0, mean = inv_coeff, sd = 1))^(1-pred_label)
    lower_bound <- (pnorm(q = 0, mean = inv_coeff, sd = 1))*pred_label
    prob_l <- (U*(upper_bound - lower_bound))+lower_bound
    prob_l[prob_l == 1] <- (1 - epsilon)
    prob_l[prob_l == 0] <- epsilon
    lat_var <- qnorm(prob_l, mean = inv_coeff, sd = 1)
    
# ========================== Generating mix_var ================================    
    wav_lat_var <- wd(lat_var, filter.number = number_vm, family =family_choice)
    lat_var_coeff <- c(1:length(lat_var))
    lat_var_coeff[1] <- wav_lat_var$C[2^(log2(length(lat_var))+1)-1]
    lat_var_coeff[2:length(lat_var)] <-wav_lat_var$D
    
    r <- 0
    for (jj in nlevelsWT(wav_lat_var):1){
      q <- 2^(jj-1)
      ninho1 <- sum(mix_var[(2+r):(r+q+1)])
      ninho0 <- q - ninho1
      parameters_pi <- rbeta(1, shape1 = (zeta+ninho1),
                             shape2 = (rho+ninho0))
      
      parameters_tau <- rgamma(1,
                               shape = (kappa+(q/2)),
                               rate = (xi+
                                         sum(coef_wavelet[(2+r):(r+q+1)]^2)/2))
      
      pi_mix_wavelet[(2+r):(r+q+1)] <- rep_len(parameters_pi,q)
      
      tau2_mix_gaussian[(2+r):(r+q+1)] <- rep_len(parameters_tau,q)
      
      r <- r + q
    }
    
    tau2_mix_gaussian <- 1/tau2_mix_gaussian
    lambda_mix_gaussian <- (tau2_mix_gaussian)/(tau2_mix_gaussian+1)
    
    pi_post_mix_wavelet <- wpost_normal(w = pi_mix_wavelet,
                                        xx = lat_var_coeff,
                                        tau2 = tau2_mix_gaussian)
    
    mix_var <- rbinom(n = simulationsize, size = 1, prob = pi_post_mix_wavelet)
    
# ============================ Generating theta ================================  
    coef_wavelet <- mix_var*rnorm(n = simulationsize,
                                  mean = lambda_mix_gaussian*
                                    lat_var_coeff,
                                  sd = sqrt(lambda_mix_gaussian))
    
    
    coef_wavelet[1] <- accessC.wd(wav_lat_var, level = 0)
# =========================== Generating inv_coeff =============================
    inv_coeff <- inverse_coeff(coef_wavelet, vm = number_vm, 
                               fc = family_choice)
    
# ========================== Generating alpha ==================================
    alpha[ii,] <- pnorm(inv_coeff, 0, 1)
    
    ### heavy calculations
    if( ii %% gc_after_a_number == 0){
      gc()
    }
  }
  
# =========================== Burn-in and lag ==================================
  #Creating the matrices with the burn-in and lag
  
  matrix_alpha <- alpha[N,]
  matrix_mu <-  mu[N,]
  matrix_phi <- phi[N,]
  
  
  estimated_alpha_median[pp,] <- apply(matrix_alpha,2,median)
  estimated_alpha_mean[pp,] <- apply(matrix_alpha,2,mean)
  estimated_mu_median[pp,] <-   apply(matrix_mu,2,median)
  estimated_tau_median[pp,] <-  apply(matrix_phi,2,median)
}

estimativa_mediana <- apply(estimated_alpha_median,2,mean)
estimativa_media <- apply(estimated_alpha_mean,2,mean)
estimate_mu <- apply(estimated_mu_median,2,mean)
estimate_tau <- apply(estimated_tau_median,2,mean)


#Final plot
plot(seq(1,simulationsize,1),alpha_behavior,type="l",col="red", lwd = 2, 
     xlab="t", ylab=expression(alpha[t]))
lines(seq(1,simulationsize,1), estimativa_mediana, lty=2, lwd = 2)
lines(seq(1,simulationsize,1), estimativa_media, lty=3,col = "blue",lwd = 3)

#Component parameters
#mu1
HDI_mu1 <-  ComputeHDI(estimated_mu_median[,1], credible.region = 0.95)
HDI_tau1 <- ComputeHDI(estimated_tau_median[,1], credible.region = 0.95)
HDI_mu2 <-  ComputeHDI(estimated_mu_median[,2], credible.region = 0.95)
HDI_tau2 <- ComputeHDI(estimated_tau_median[,2], credible.region = 0.95)

#Component parameters
#mu1
round(estimate_mu[1], digits=2)
round(HDI_mu1, digits=2)

#tau1
round(estimate_tau[1], digits=2)
round(HDI_tau1, digits=2)

#mu2
round(estimate_mu[2], digits=2)
round(HDI_mu2, digits=2)

#tau2
round(estimate_tau[2], digits=2)
round(HDI_tau2, digits=2)

