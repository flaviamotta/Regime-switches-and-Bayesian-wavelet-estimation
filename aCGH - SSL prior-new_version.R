# ==============================================================================
# ================================= Packages ===================================
# ==============================================================================
#aCGH - SSL prior

library("wavethresh")
library("EbayesThresh")
library('waveslim')
library('changepoint')
library('ggplot2')
library('latex2exp')
library('coda')
library('ggtext')
library('bfw')

# ==============================================================================
# ============================= Preparing the data =============================
# ==============================================================================


data(Lai2005fig4)#data(Lai2005fig3)
y <- Lai2005fig4$GBM29#y <- Lai2005fig3$GBM31
data_size <- length(y)
aux_data_size <- 2^ceiling(log2(data_size ))
wavelet_adapted_size <- aux_data_size - data_size
serie <- seq_len(data_size)/data_size

# ============================ Global variables ================================

#Gibbs sampler: Setting down the size of our chain, of our burn-in and our lag

burn <- 1000
lags <- 50
nchain <- 1000
BB <- burn + lags*nchain
n.verbose <- 1000
verbose <- TRUE

#Family choice
family_choice <- "DaubExPhase" #"DaubLeAsymm", "DaubExPhase", "Coiflets"

#Number of vanishing moments
number_vm <- 10


#Function to transform the wavelet coefficients back to the data domain
inverse_coeff <- function(coeff, vm, fc){
  zero <- rep(0, length(coeff))
  zero_wd <- wd(zero, filter.number = vm, family = fc)
  zero_wd$C[2^(log2(length(coeff))+1)-1] <- coeff[1]
  zero_wd$D <- coeff[2:length(coeff)]
  return(wr(zero_wd))    
}

# Function to calculate the posterior weight for non-null effect
wpost_laplace <- function(w, x, s = 1, a){
  wpostt <- 1 - (1 - w)/(1 + w * beta.laplace(x, s, a))
  wpostt[1] <- 1
  return(wpostt)
}

# ==============================================================================
# ============================== Gibbs Sampler =================================
# ==============================================================================

mu <- matrix(nrow=BB,ncol = 2)
phi <- matrix(nrow=BB,ncol = 2)
coef_wavelet <- matrix(nrow=BB,ncol=aux_data_size)
alpha <- matrix(nrow=BB,ncol=data_size) 
inv_coeff <- matrix(nrow=BB,ncol=aux_data_size)

prior_mean1 <- quantile(y, probs = 0.25)
prior_mean2 <- quantile(y, probs = 0.75)
prior_variance1 <- 10*var(y)
prior_variance2 <- 10*var(y)
prior_alpha_gamma <- 0.01
prior_beta_gamma <- 0.01

#Defining our start values
mu[1,]<-c(prior_mean1,prior_mean2)
phi[1,]<- c(1/prior_variance1, 1/prior_variance2)
coef_wavelet[1,] <- vector("double", aux_data_size) 
alpha[1,] <- rep_len(1, data_size)  
prob0start <- (1-alpha[1,])*dnorm(y, mean = mu[1,1], 
                                  sd = 1/sqrt(phi[1,1]))
prob1start <- alpha[1,]*dnorm(y, mean = mu[1,2], 
                              sd = 1/sqrt(phi[1,2]))
probvalstart <- rep_len(1, data_size)

pred_label <- rbinom(data_size, 1, prob = probvalstart)
inv_coeff[1,] <- inverse_coeff(coef_wavelet[1,], vm = number_vm,
                               fc = family_choice) 
lat_var <- rnorm(aux_data_size, mean = inv_coeff[1,], 1)
mix_var <- rbinom(n = aux_data_size, size = 1, prob = 0.1)
ninho1 <- sum(mix_var)
ninho0 <- data_size - ninho1
pi_mix_wavelet <- rep_len(1, data_size)
a_laplace <- rep_len(1, data_size)
epsilon <- 1e-10  # Small constant to prevent division by zero

zeta <- 1
rho <- 1
kappa <- 1
xi <- 100


# ================================== Gibbs =====================================

for(ii in 2:BB){
  #Creating some auxiliary variables for calculating the conditional posterior 
  n1 <- sum(pred_label)
  n0 <- data_size - n1
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
  #For dealing with label switching, we stablish that the pairs (\mu_k, \phi_k) 
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
  
  pred_label <- rbinom(n = data_size, size = 1, prob = probval)
  
# ========================= Generating lat_var =================================
  pred_label_aux <- c(pred_label,rev(pred_label)[1:wavelet_adapted_size])
  U <- runif(n = aux_data_size)
  upper_bound <-(pnorm(q = 0,mean =inv_coeff[ii-1,], sd = 1))^(1-pred_label_aux)
  lower_bound <- (pnorm(q = 0, mean = inv_coeff[ii-1,], sd = 1))*pred_label_aux
  prob_l <- (U*(upper_bound - lower_bound))+lower_bound
  prob_l[prob_l == 1] <- (1 - epsilon)
  prob_l[prob_l == 0] <- epsilon
  lat_var <- qnorm(prob_l, mean = inv_coeff[ii-1,], sd = 1)
  
# ========================== Generating mix_var ================================    
  wav_lat_var <- wd(lat_var, filter.number = number_vm, family = family_choice)
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
    parameters_a <- rgamma(1,
                           shape = (kappa+q),
                           rate = (xi+
                                     sum(abs(coef_wavelet[ii-1,(2+r):(r+q+1)])))
    )
    
    pi_mix_wavelet[(2+r):(r+q+1)] <- rep_len(parameters_pi,q)
    
    a_laplace[(2+r):(r+q+1)] <- rep_len(parameters_a,q)
    
    r <- r + q
  }
  
  pi_post_mix_wavelet <- wpost_laplace(w = pi_mix_wavelet,
                                       x = lat_var_coeff,
                                       a = a_laplace)
  
  mix_var <- rbinom(n = aux_data_size, size = 1, prob = pi_post_mix_wavelet)
  
# ============================ Generating theta ================================  
  prob_tn_nu <- exp(-a_laplace*lat_var_coeff)*pnorm((lat_var_coeff-a_laplace))
  
  prob_tn_de <- prob_tn_nu +
    exp(a_laplace*lat_var_coeff)*(1-pnorm((lat_var_coeff+a_laplace)))
  
  prob_tn <- prob_tn_nu / (prob_tn_de + epsilon)
  
  pre_tn <- rbinom(n = aux_data_size, size = 1, prob = prob_tn)
  
  
  mean_tn <- (lat_var_coeff-a_laplace)*pre_tn +
    (lat_var_coeff+a_laplace)*(1-pre_tn)
  
  
  upper_bound_tn <- (pnorm(q = 0, mean = mean_tn, sd = 1))^(1-pre_tn)
  lower_bound_tn <- (pnorm(q = 0, mean = mean_tn, sd = 1))*pre_tn
  
  prob_l_tn <- lower_bound_tn + 
    runif(n = aux_data_size) * (upper_bound_tn - lower_bound_tn)
  prob_l_tn[prob_l_tn == 1] <- (1 - epsilon)
  prob_l_tn[prob_l_tn == 0] <- epsilon
  coef_wavelet[ii,] <- mix_var*qnorm(prob_l_tn, mean = mean_tn, sd = 1)
  
  coef_wavelet[ii,1] <- accessC.wd(wav_lat_var, level = 0)  
  
# =========================== Generating inv_coeff =============================
  inv_coeff[ii,] <- inverse_coeff(coef_wavelet[ii,], vm = number_vm, 
                                  fc = family_choice)
  
# ========================== Generating alpha ==================================
  alpha[ii,] <- pnorm(inv_coeff[ii,], 0, 1)[1:data_size]
  
  if(!ii%%n.verbose & verbose)
    print(paste0(ii, " iterations of ", BB, "."))
}




# =========================== Burn-in and lag ==================================
#Creating the matrices with the burn-in and lag

N <- seq(burn+1, BB, lags)

matrix_mu <- mu[N,]
matrix_phi <- phi[N,]
matrix_alpha <- alpha[N,]

#estimated values by median
estimated_mu_median <- apply(matrix_mu,2,median)
estimated_phi_median <- apply(matrix_phi,2,median)
estimated_alpha_median <- apply(matrix_alpha,2,median)

#estimated values by mean
estimated_mu_mean <- apply(matrix_mu,2,mean)
estimated_phi_mean <- apply(matrix_phi,2,mean)
estimated_alpha_mean <- apply(matrix_alpha,2,mean)


# ================================== Plots =====================================


plot(estimated_alpha_median,type="l",col="black", lwd = 1, 
     ylab=expression(alpha[t]))
lines(seq(1,data_size,1), estimated_alpha_mean, lty=3, lwd = 3,col = "blue")

#Trace mu1
plot.ts(matrix_mu[,1], main = "", xlab="", 
        ylab = expression(mu[1]))
abline(h = estimated_mu_median[1], lty=2, col = "orange", lwd = 3)
abline(h = estimated_mu_mean[1], lty=3, col = "blue", lwd = 3)

#Correlogram mu1
acf(matrix_mu[,1], main = expression(mu[1]))

#Density mu1
plot(density(matrix_mu[,1]), main = "", col="black", lwd = 2, xlab="", 
     ylab = expression(paste("Density ",mu[1])))
abline(v = estimated_mu_median[1], lty=2, col = "orange", lwd = 3)
abline(v = estimated_mu_mean[1], lty=3, col = "blue", lwd = 3)


# ==============================================================================

#Trace mu2
plot.ts(matrix_mu[,2], main = "", xlab="", 
        ylab = expression(mu[2]))
abline(h = estimated_mu_median[2], lty=2, col = "orange", lwd = 3)
abline(h = estimated_mu_mean[2], lty=3, col = "blue", lwd = 3)


#Correlogram mu2
acf(matrix_mu[,2], main = expression(mu[2]))


#Density mu2
plot(density(matrix_mu[,2]), main = "", col="black", lwd = 2, xlab="", 
     ylab = expression(paste("Density ",mu[2])))
abline(v = estimated_mu_median[2], lty=2, col = "orange", lwd = 3)
abline(v = estimated_mu_mean[2], lty=3, col = "blue", lwd = 3)

# ==============================================================================

#Trace phi1
plot.ts(matrix_phi[,1], main = "", xlab="", 
        ylab = expression(phi[1]))
abline(h = estimated_phi_median[1], lty=2, col = "orange", lwd = 3)
abline(h = estimated_phi_mean[1], lty=3, col = "blue", lwd = 3)

#Correlogram phi1
acf(matrix_phi[,1], main = expression(phi[1]))

#Density phi1
plot(density(matrix_phi[,1]), main = "", col="black", lwd = 2, xlab="", 
     ylab = expression(paste("Density ",phi[1])))
abline(v = estimated_phi_median[1], lty=2, col = "orange", lwd = 3)
abline(v = estimated_phi_mean[1], lty=3, col = "blue", lwd = 3)

# ==============================================================================

#Trace phi2
plot.ts(matrix_phi[,2], main = "", xlab="", 
        ylab = expression(phi[2]))
abline(h = estimated_phi_median[2], lty=2, col = "orange", lwd = 3)
abline(h = estimated_phi_mean[2], lty=3, col = "blue", lwd = 3)


#Correlogram phi2
acf(matrix_phi[,2], main = expression(phi[2]))


#Density phi2
plot(density(matrix_phi[,2]), main = "", col="black", lwd = 2, xlab="", 
     ylab = expression(paste("Density ",phi[2])))
abline(v = estimated_phi_median[2], lty=2, col = "orange", lwd = 3)
abline(v = estimated_phi_mean[2], lty=3, col = "blue", lwd = 3)

# ==============================================================================

#Alpha

minHPD <- vector("double", data_size)
maxHPD <- vector("double", data_size)
for (ii in 1:data_size) {
  hpd_alpha <- ComputeHDI(matrix_alpha[,ii], credible.region = 0.95)
  minHPD[ii] <- hpd_alpha[1]
  maxHPD[ii] <- hpd_alpha[2]
}

alpha_post <- data.frame(minHPD, maxHPD, estimated_alpha_median,
                         serie, y)

ggplot(alpha_post, aes(x = serie))+
  geom_text(aes(x = 0.5, y = 0.5,
                label = ""),
            stat = "unique") +
  geom_text(aes(x = 0.5, y = 1,
                label = ""),
            stat = "unique") +
  geom_ribbon(aes(ymin = minHPD, ymax = maxHPD), alpha = .75,
              fill = "grey", color = "transparent")+
  geom_line(aes(y = estimated_alpha_median), color = "black", size = 1)+
  labs(title="",x="t", y = TeX("$\\hat{\\alpha}_t$")) +
  theme_classic()


ggplot(alpha_post, aes(x = serie))+
  geom_text(aes(x = 0.5, y = 0.5,
                label = ""),
            stat = "unique") +
  geom_text(aes(x = 0.5, y = 1,
                label = ""),
            stat = "unique") +
  geom_line(aes(y = estimated_alpha_median), color = "black", size = 1)+
  labs(title="",x="t", y = TeX("$\\hat{\\alpha}_t$")) +
  theme_classic()

