#' Equilibrium initialisation list creation: P. vivax
#'
#' The function \code{vivax_equilibrium_init_create} creates an equilibrium initialisation state to be
#' used within later model runs
#'
#' @param age Vector of age brackets.
#' @param het_brackets Integer number of biting heteogenity compartments.
#' @param country String for country of interest. If NULL the seasonal parameters
#' will attempt to be loaded using just the admin unit, however if there is ambiguity
#' in the admin unit an error will be thrown. If both NULL then no seasonality is
#' assumed. Default = NULL.
#' @param admin_unit String for admin unit with country for loading seasonal
#' parameters. If country is NULL, the admin unit will attempt to be located,however
#' if there is ambiguity in the admin unit an error will be thrown. If both country
#' and admin_unit are NULL then no seasonality is assumed. Default = NULL.
#' @param ft Numeric for the frequency of people seeking treatment.
#' @param EIR Numeric for desired annual EIR.
#' @param K_max Maximum number of hypnozoite batches. Has to be an integer between
#' 1 and 30. Note a model with stratification into larger number of hypnozoite batches
#' (e.g. K_max > 5) may have to be run on the cluster. For an EIR of 1 or less,
#' K_max = 2 is sufficient.
#' @param model_param_list List of epidemiological parameters created by
#'
#' @importFrom stringi stri_trans_general
#' @importFrom statmod gauss.quad.prob
#'
#'
#' @export

vivax_equilibrium_init_create_combined <- function(age, ft,
                                          EIR, p, K_max = 30, 
                                          
                                          MW_age_rates_prop = NULL, 
                                          use_mid_ages = NULL,
                                          divide_omega = NULL,
                                          do_not_integrate = NULL,
                                          malariasimulationoutput = F)
{
  
  ## Check Parameters
  if(!is.numeric(age)) stop("age provided is not numeric")
  if(!is.numeric(p$n_heterogeneity_groups)) stop("het_brackets provided is not numeric")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  if(!is.numeric(EIR)) stop("EIR provided is not numeric")
  if(!is.numeric(K_max)) stop("K_max is not provided or not numeric")
  if(!(K_max %in% c(1:30))) stop("K_max is not an integer between 1 and 30")
  
  p$dp <- 10
  # browser()
  ## Handle parameters
  ## Population demographics
  age <- age * 365
  max_age <- 100*365
  p$eta <- 1/p$average_age
  na <- as.integer(length(age)-1)  # number of age groups
  nh <- as.integer(p$n_heterogeneity_groups) # number of heterogeneity groups
  nk <- as.integer(K_max+1)      # number of hypnozoite batch groups
  h <- statmod::gauss.quad.prob(nh, dist = "normal")
  
  ## Hypnozoite inputs
  K_max_switch <- c(rep(1,nk-1),0)        # switch to set some transitions to 0 when k=K_max
  K_max_switch_on <- c(rep(0,nk-1),1)     # switch to only apply some transitions when k=K_max
  K0_switch <- c(0, rep(1,nk-1))          # switch to set some transitions to 0 when k=0
  kk_vec <- 0:K_max                       # vector of hypnozoite batches
  
  age_rate <- age_width <- age_mid_point <- den <- c()
  
  if(is.null(MW_age_rates_prop)){
    
    for (i in 1:(na-1))
    {
      age_width[i] <- age[i+1] - age[i]
      age_rate[i] <- 1/(age[i + 1] - age[i])  # vector of rates at which people leave each age group (1/age group width)
      age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group
      
    }
    age_rate[na] <- 0
    # age_mid_point[na] <- 0.5 * (age[na] + max_age)
    
    age_mid_point[na] <- age[na]
    
    den <- 1/(1 + age_rate[1]/p$eta)
    for (i in 1:(na-1))
    {
      den[i+1] <- age_rate[i] * den[i]/(age_rate[i+1] + p$eta)  # proportion in each age group
    }
    
    
    
  } else {
    
    den <- rep(NA, na)
    den[1:(na-1)]  <-  exp( - age[1:(na-1)]*p$eta) - exp( - age[2:na]*p$eta )
    den[na]        <- 1 - sum( den[1:(na-1)] )
    
    ## Calculate rate of aging, 1/duration in age group
    ## r is equivalent to r in MW eq
    r <- rep(NA, na)
    for(i in 1:(na-1)){
      r[i] <- ( 1 - sum(den[1:i]) )/( den[i]/p$eta )
      age_width[i] <- age[i+1] - age[i]
      age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group
    }
    
    r[na] = 0
    
    age_rate <- rep(NA, na)
    age_rate[1:(na-1)] <- r[1:(na-1)]*(den[1:(na-1)]/den[2:na])
    age_rate[na] <- 0
    # age_mid_point[na] <- 0.5 * (age[na] + max_age)
    age_mid_point[na] <- age[na]
    
    
  }
  
  ## force of infection
  foi_age <- c()
  if(is.null(use_mid_ages)){
    for (i in 1:na)
    {
      foi_age[i] <- 1 - (p$rho * exp(-age[i]/p$a0))  #force of infection for each age group
    }
  } else {
    for (i in 1:na)
    {
      foi_age[i] <- 1 - (p$rho * exp(-age_mid_point[i]/p$a0))  #force of infection for each age group
    }
  }
  ## NOTE: in vivax package foi_age uses age_mid_point instead of age
  fden <- foi_age * den
  omega <- sum(fden)  #normalising constant
  
  if(!is.null(divide_omega)){
    foi_age <- foi_age/omega
  }
  
  
  
  ## heterogeneity
  het_x <- h$nodes
  het_wt <- h$weights
  den_het <- outer(den, het_wt)
  rel_foi <- exp(-p$sigma_squared/2 + sqrt(p$sigma_squared) * het_x)/sum(het_wt * exp(-p$sigma_squared/2 + sqrt(p$sigma_squared) * het_x))
  # rel_foi <- rel_foi/omega
  ## EIR
  EIRY_eq <- EIR  # initial annual EIR
  EIRd_eq <- EIRY_eq/365
  EIR_eq <- outer(foi_age, rel_foi) * EIRd_eq
  
  # maternal immunity begins at a level proportional to the clinical
  # immunity of a 20 year old, this code finds that level
  age20i <- rep(0, na)
  for (i in 2:na)
  {
    age20i[i] <- ifelse(age[i] >= (20 * 365) & age[i - 1] < (20 * 365),
                        i, age20i[i - 1])
  }
  age20u <- as.integer(age20i[na])
  age20l <- as.integer(age20u - 1)
  age_20_factor <- (20 * 365 - age[age20l] - 0.5 * age_width[age20l]) *
    2/(age_width[age20l] + age_width[age20u])
  
  
  # Find equilibrium FOI
  FOI_eq <- matrix(NA, nrow=na, ncol=nh)
  
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      FOI_eq[i, j] <- EIR_eq[i, j] * p$b
    }
  }
  
  
  if(is.null(do_not_integrate)){
  
  ## Equilibrium state values for humans
  
  # Isolate demographic component in common to all human equations:
  gamma <- p$eta + c(age_rate[1:(na - 1)], 0)   # leaving current age group (death + aging out)
  delta <- c(p$eta, age_rate[1:(na - 1)])       # entering current age group (aging in)
  # At age 0, rate of aging in = mortality rate, to replace population
  
  # All equilibrium state values are calculated using matrix inversion (solve function)
  # Notation: A * X + b = 0
  # where A is a matrix of transitions in differential equations
  
  ## Equilibrium distribution of hypnozoite batches
  Hyp_eq <- array(NA, dim=c(na,nh,nk))
  
  b_Hyp <- rep(list(array(0,dim=c(1,nk,nh))), na)
  
  # Create matrix A for each age and heterogeneity group
  A_Hyp <- rep(list(array(0,dim=c(nk,nk,nh))), na)
  
  for(j in 1:nh) {
    for(i in 1:na) {
      for (kk in 1:(nk-1)) {
        A_Hyp[[i]][kk,kk,j] <- -FOI_eq[i,j]-kk_vec[kk]*p$gammal-gamma[i]
        A_Hyp[[i]][kk,kk+1,j] <- kk_vec[kk+1]*p$gammal
        A_Hyp[[i]][kk+1,kk,j] <- FOI_eq[i,j]
      }
      A_Hyp[[i]][nk,nk,j] <- -K_max*p$gammal-gamma[i]
      
    }
  }
  
  # Calculate equilibrium in each hypnozoite state at age 0:
  for(j in 1:nh) {
    b_Hyp[[1]][,,j] <-  c(-p$eta*het_wt[j],rep(0,nk-1))
    
    # Solve matrices:
    Hyp_eq[1,j,] <- solve(A_Hyp[[1]][,,j], b_Hyp[[1]][,,j])
  }
  
  # For all other age groups:
  for(j in 1:nh) {
    for(i in 2:na) {
      b_Hyp[[i]][,,j] <-  c(-delta[i]*Hyp_eq[i-1,j,])
      
      # Solve matrices:
      Hyp_eq[i,j,] <- solve(A_Hyp[[i]][,,j], b_Hyp[[i]][,,j])
      
    }
  }
  
  # Calculate proportion in each hypnozoite state for each age and heterogeneity group
  hyp_wt <- array(NA, dim=c(na,nh,nk))
  for (j in 1:nh) {
    for (i in 1:na) {
      for (k in 1:nk) {
        hyp_wt[i,j,k] <- Hyp_eq[i,j,k]/sum(Hyp_eq[i,j,])
      }
    }
  }
  # browser()
  # hyp_wt
  # sum(colSums(Hyp_eq, dims = 2)*0:10)
  
  
  ## Immunity states
  AP_eq <- array(NA, dim=c(na,nh,nk))
  AC_eq <- array(NA, dim=c(na,nh,nk))
  AP_MAT_init_eq <- vector(length = nh, mode = "numeric")
  AP_MAT_eq <- matrix(0, na, nh)
  AC_MAT_init_eq <- vector(length = nh, mode = "numeric")
  AC_MAT_eq <- matrix(0, na, nh)
  
  b_AP <- list()
  b_AC <- list()
  
  # Create matrix A for each age and heterogeneity group for anti-parasite (AP) and anti-clinical immunity (AC)
  A_AP <- rep(list(array(0,dim=c(nk,nk,nh))), na)
  A_AC <- rep(list(array(0,dim=c(nk,nk,nh))), na)
  
  for(j in 1:nh) {
    for(i in 1:na) {
      for (kk in 1:(nk-1)) {
        A_AP[[i]][kk,kk,j] <- -FOI_eq[i,j]-1/p$rid-(kk-1)*p$gammal-gamma[i]
        A_AP[[i]][kk,kk+1,j] <- p$gammal*kk * Hyp_eq[i,j,kk+1]/Hyp_eq[i,j,kk]
        A_AP[[i]][kk+1,kk,j] <- FOI_eq[i,j] * Hyp_eq[i,j,kk]/Hyp_eq[i,j,kk+1]
        
        A_AC[[i]][kk,kk,j] <- -FOI_eq[i,j]-1/p$rc-(kk-1)*p$gammal-gamma[i]
        A_AC[[i]][kk,kk+1,j] <- p$gammal*kk * Hyp_eq[i,j,kk+1]/Hyp_eq[i,j,kk]
        A_AC[[i]][kk+1,kk,j] <- FOI_eq[i,j] * Hyp_eq[i,j,kk]/Hyp_eq[i,j,kk+1]
      }
      A_AP[[i]][nk,nk,j] <- -1/p$rid-K_max*p$gammal-gamma[i]
      A_AC[[i]][nk,nk,j] <- -1/p$rc-K_max*p$gammal-gamma[i]
      
    }
  }
  
  # Calculate equilibrium in each AP and AC hypnozoite state at age 0:
  b_AP <- rep(list(array(0,dim=c(1,nk,nh))), na)
  b_AC <- rep(list(array(0,dim=c(1,nk,nh))), na)
  
  for(j in 1:nh) {
    for(kk in 2:nk) {
      b_AP[[1]][,kk,j] <- -(FOI_eq[1,j]*Hyp_eq[1,j,kk-1]/Hyp_eq[1,j,kk] + kk_vec[kk]*p$f)/
        ((FOI_eq[1,j]*Hyp_eq[1,j,kk-1]/Hyp_eq[1,j,kk] + kk_vec[kk]*p$f)*p$ud+1)
      
      b_AC[[1]][,kk,j] <- -(FOI_eq[1,j]*Hyp_eq[1,j,kk-1]/Hyp_eq[1,j,kk] + kk_vec[kk]*p$f)/
        ((FOI_eq[1,j]*Hyp_eq[1,j,kk-1]/Hyp_eq[1,j,kk] + kk_vec[kk]*p$f)*p$uc+1)
      
    }
    
    # Solve matrices:
    AP_eq[1,j,] <- solve(A_AP[[1]][,,j], b_AP[[1]][,,j])
    AC_eq[1,j,] <- solve(A_AC[[1]][,,j], b_AC[[1]][,,j])
  }
  
  # For all other age groups:
  for(j in 1:nh) {
    for(i in 2:na) {
      
      b_AP[[i]][,1,j] <- -delta[i]*AP_eq[i-1,j,1]*Hyp_eq[i-1,j,1]/Hyp_eq[i,j,1]
      b_AC[[i]][,1,j] <- -delta[i]*AC_eq[i-1,j,1]*Hyp_eq[i-1,j,1]/Hyp_eq[i,j,1]
      
      for(kk in 2:nk) {
        b_AP[[i]][,kk,j] <- -(FOI_eq[i,j]*Hyp_eq[i,j,kk-1]/Hyp_eq[i,j,kk] + kk_vec[kk]*p$f)/
          ((FOI_eq[i,j]*Hyp_eq[i,j,kk-1]/Hyp_eq[i,j,kk] + kk_vec[kk]*p$f)*p$ud+1)-
          delta[i]*AP_eq[i-1,j,kk]*Hyp_eq[i-1,j,kk]/Hyp_eq[i,j,kk]
        
        b_AC[[i]][,kk,j] <- -(FOI_eq[i,j]*Hyp_eq[i,j,kk-1]/Hyp_eq[i,j,kk] + kk_vec[kk]*p$f)/
          ((FOI_eq[i,j]*Hyp_eq[i,j,kk-1]/Hyp_eq[i,j,kk] + kk_vec[kk]*p$f)*p$uc+1)-
          delta[i]*AC_eq[i-1,j,kk]*Hyp_eq[i-1,j,kk]/Hyp_eq[i,j,kk]
      }
      
      # Solve matrices:
      AP_eq[i,j,] <- solve(A_AP[[i]][,,j], b_AP[[i]][,,j])
      AC_eq[i,j,] <- solve(A_AC[[i]][,,j], b_AC[[i]][,,j])
      
    }
  }
  
  # Maternal immunity needs to be calculated after because it references AP and AC
  
  # Immunity in 20 year old woman is averaged over hypnozoite states
  # (weighted by hypnozoite distribution)
  
  AP_eq_mean <- apply(AP_eq*hyp_wt, c(1,2), sum)
  AC_eq_mean <- apply(AC_eq*hyp_wt, c(1,2), sum)
  
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      AP_MAT_init_eq[j] <- p$pcm * (AP_eq_mean[age20l, j] + age_20_factor *
                                      (AP_eq_mean[age20u, j] - AP_eq_mean[age20l, j]))
      AC_MAT_init_eq[j] <- p$pcm * (AC_eq_mean[age20l, j] + age_20_factor *
                                      (AC_eq_mean[age20u, j] - AC_eq_mean[age20l, j]))
      
      AP_MAT_eq[i, j] <- AP_MAT_init_eq[j]*exp(-age[i]/p$rm)
      AC_MAT_eq[i, j] <- AC_MAT_init_eq[j]*exp(-age[i]/p$rm)
      # Note IBM uses age_mid_point
      
    }
  }
  
  
  # Turn AP_MAT_eq and AC_MAT_eq into an array (with same values for each hypnozoite state)
  # for calculation of phi_LM, phi_D and dPCR
  AP_MAT_eq <- array(AP_MAT_eq, dim=c(na,nh,(K_max+1)))
  AC_MAT_eq <- array(AC_MAT_eq, dim=c(na,nh,(K_max+1)))
  
  # Calculate probabilities
  phi_LM_eq <- p$philm_min + (p$philm_max-p$philm_min) *
    1/(1+((AP_eq+AP_MAT_eq)/p$alm50)^p$klm)
  phi_D_eq <- p$phi0*p$phi1 + (p$phi0 - p$phi0*p$phi1) *
    1/(1+((AC_eq+AC_MAT_eq)/p$ic0)^p$kc)
  dPCR_eq <- p$dpcr_min + (p$dpcr_max-p$dpcr_min) *
    1/(1+((AP_eq+AP_MAT_eq)/p$apcr50)^p$kpcr)
  
  ## Human states (stored in Z)
  n <- 6*nk       # matrix dimensions = epidemiological compartments x hypnozoite groups
  Z <- array(NA, dim=c(na,nh,n))
  b <- list()
  
  # Define indices:
  iS <- 1
  iI_PCR <- 2
  iI_LM <- 3
  iI_D <- 4
  iT <- 5
  iP <- 6
  
  # Create matrix A for each age and heterogeneity group
  A <- array(0,dim=c(n,n,nh))
  A <- rep(list(A), na)
  
  for(j in 1:nh) {
    for(i in 1:na) {
      
      # Fill in non-zero elements:
      
      # Transitions for k<K_max
      for (kk in 0:(K_max-1)) {
        
        # Diagonals
        A[[i]][iS+kk*6,iS+kk*6,j] <- -FOI_eq[i,j]-kk*p$f-kk*p$gammal-gamma[i]
        
        A[[i]][iI_PCR+kk*6,iI_PCR+kk*6,j] <- -FOI_eq[i,j]-kk*p$f-1/dPCR_eq[i,j,(kk+1)]-
          kk*p$gammal-gamma[i]+kk*p$f*(1-phi_LM_eq[i,j,(kk+1)])
        
        A[[i]][iI_LM+kk*6,iI_LM+kk*6,j] <- -FOI_eq[i,j]-kk*p$f-1/p$da-
          kk*p$gammal-gamma[i]+kk*p$f*(1-phi_D_eq[i,j,(kk+1)])
        
        A[[i]][iI_D+kk*6,iI_D+kk*6,j] <- -FOI_eq[i,j]-1/p$dd-kk*p$gammal-gamma[i]
        
        A[[i]][iT+kk*6,iT+kk*6,j] <- -FOI_eq[i,j]-1/p$dt-kk*p$gammal-gamma[i]
        
        A[[i]][iP+kk*6,iP+kk*6,j] <- -FOI_eq[i,j]-1/p$dp-kk*p$gammal-gamma[i]
        
        # Additions from gammal(k+1):
        A[[i]][iS+kk*6,iS+kk*6+6,j] <- (kk+1)*p$gammal
        A[[i]][iI_PCR+kk*6,iI_PCR+kk*6+6,j] <- (kk+1)*p$gammal
        A[[i]][iI_LM+kk*6,iI_LM+kk*6+6,j] <- (kk+1)*p$gammal
        A[[i]][iI_D+kk*6,iI_D+kk*6+6,j] <- (kk+1)*p$gammal
        A[[i]][iT+kk*6,iT+kk*6+6,j] <- (kk+1)*p$gammal
        A[[i]][iP+kk*6,iP+kk*6+6,j] <- (kk+1)*p$gammal
        
        # Other transitions:
        A[[i]][iI_PCR+kk*6,iS+kk*6,j] <- kk*p$f*(1-phi_LM_eq[i,j,(kk+1)])
        A[[i]][iI_LM+kk*6,iS+kk*6,j] <- kk*p$f*(1-phi_D_eq[i,j,(kk+1)])*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iI_LM+kk*6,iI_PCR+kk*6,j] <- kk*p$f*(1-phi_D_eq[i,j,(kk+1)])*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iI_D+kk*6,iS+kk*6,j] <- kk*p$f*phi_D_eq[i,j,(kk+1)]*(1-ft)*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iI_D+kk*6,iI_PCR+kk*6,j] <- kk*p$f*phi_D_eq[i,j,(kk+1)]*(1-ft)*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iI_D+kk*6,iI_LM+kk*6,j] <- kk*p$f*phi_D_eq[i,j,(kk+1)]*(1-ft)
        A[[i]][iT+kk*6,iS+kk*6,j] <- kk*p$f*phi_D_eq[i,j,(kk+1)]*ft*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iT+kk*6,iI_PCR+kk*6,j] <- kk*p$f*phi_D_eq[i,j,(kk+1)]*ft*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iT+kk*6,iI_LM+kk*6,j] <- kk*p$f*phi_D_eq[i,j,(kk+1)]*ft
        
      }
      
      # Diagonals for K_max:
      A[[i]][iS+K_max*6,iS+K_max*6,j] <- -FOI_eq[i,j]-K_max*p$f-K_max*p$gammal-gamma[i]
      A[[i]][iI_PCR+K_max*6,iI_PCR+K_max*6,j] <- -FOI_eq[i,j] + FOI_eq[i,j]*(1-phi_LM_eq[i,j,(K_max+1)]) -
        K_max*p$f-1/dPCR_eq[i,j,(K_max+1)]-
        K_max*p$gammal-gamma[i]+K_max*p$f*(1-phi_LM_eq[i,j,(K_max+1)])
      A[[i]][iI_LM+K_max*6,iI_LM+K_max*6,j] <- -FOI_eq[i,j] + FOI_eq[i,j]*(1-phi_D_eq[i,j,(K_max+1)])-
        K_max*p$f-1/p$da-
        K_max*p$gammal-gamma[i]+K_max*p$f*(1-phi_D_eq[i,j,(K_max+1)])
      A[[i]][iI_D+K_max*6,iI_D+K_max*6,j] <- -1/p$dd-K_max*p$gammal-gamma[i]
      A[[i]][iT+K_max*6,iT+K_max*6,j] <- -1/p$dt-K_max*p$gammal-gamma[i]
      A[[i]][iP+K_max*6,iP+K_max*6,j] <- -1/p$dp-K_max*p$gammal-gamma[i]
      
      # Other transitions for K_max:
      A[[i]][iI_PCR+K_max*6,iS+K_max*6,j] <- FOI_eq[i,j]*(1-phi_LM_eq[i,j,(K_max+1)])+
        K_max*p$f*(1-phi_LM_eq[i,j,(K_max+1)])
      A[[i]][iI_LM+K_max*6,iS+K_max*6,j] <- FOI_eq[i,j]*phi_LM_eq[i,j,(K_max+1)]*(1-phi_D_eq[i,j,(K_max+1)])+
        K_max*p$f*(1-phi_D_eq[i,j,(K_max+1)])*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iI_LM+K_max*6,iI_PCR+K_max*6,j] <- FOI_eq[i,j]*phi_LM_eq[i,j,(K_max+1)]*(1-phi_D_eq[i,j,(K_max+1)])+
        K_max*p$f*(1-phi_D_eq[i,j,(K_max+1)])*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iI_D+K_max*6,iS+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*(1-ft)*phi_LM_eq[i,j,(K_max+1)] +
        K_max*p$f*phi_D_eq[i,j,(K_max+1)]*(1-ft)*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iI_D+K_max*6,iI_PCR+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*(1-ft)*phi_LM_eq[i,j,(K_max+1)] +
        K_max*p$f*phi_D_eq[i,j,(K_max+1)]*(1-ft)*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iI_D+K_max*6,iI_LM+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*(1-ft)+
        K_max*p$f*phi_D_eq[i,j,(K_max+1)]*(1-ft)
      A[[i]][iT+K_max*6,iS+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*ft*phi_LM_eq[i,j,(K_max+1)] +
        K_max*p$f*phi_D_eq[i,j,(K_max+1)]*ft*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iT+K_max*6,iI_PCR+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*ft*phi_LM_eq[i,j,(K_max+1)] +
        K_max*p$f*phi_D_eq[i,j,(K_max+1)]*ft*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iT+K_max*6,iI_LM+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*ft +
        K_max*p$f*phi_D_eq[i,j,(K_max+1)]*ft
      
      # Transitions for all k:
      for (kk in c(0:K_max)) {
        
        A[[i]][iS+kk*6,iI_PCR+kk*6,j] <- 1/dPCR_eq[i,j,(kk+1)]
        A[[i]][iS+kk*6,iP+kk*6,j] <- 1/p$dp
        A[[i]][iI_PCR+kk*6,iI_LM+kk*6,j] <- 1/p$da
        A[[i]][iI_LM+kk*6,iI_D+kk*6,j] <- 1/p$dd
        A[[i]][iP+kk*6,iT+kk*6,j] <- 1/p$dt
        
      }
      
      
      # Transitions for k>0 (all except S)
      for (kk in c(1:K_max)) {
        A[[i]][iI_PCR+kk*6,iS+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_LM_eq[i,j,kk])
        A[[i]][iI_PCR+kk*6,iI_PCR+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_LM_eq[i,j,kk])
        A[[i]][iI_LM+kk*6,iS+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_D_eq[i,j,kk])*phi_LM_eq[i,j,kk]
        A[[i]][iI_LM+kk*6,iI_PCR+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_D_eq[i,j,kk])*phi_LM_eq[i,j,kk]
        A[[i]][iI_LM+kk*6,iI_LM+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_D_eq[i,j,kk])
        A[[i]][iI_D+kk*6,iS+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*(1-ft)*phi_LM_eq[i,j,kk]
        A[[i]][iI_D+kk*6,iI_PCR+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*(1-ft)*phi_LM_eq[i,j,kk]
        A[[i]][iI_D+kk*6,iI_LM+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*(1-ft)
        A[[i]][iI_D+kk*6,iI_D+kk*6-6,j] <- FOI_eq[i,j]
        A[[i]][iT+kk*6,iS+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*ft*phi_LM_eq[i,j,kk]
        A[[i]][iT+kk*6,iI_PCR+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*ft*phi_LM_eq[i,j,kk]
        A[[i]][iT+kk*6,iI_LM+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*ft
        A[[i]][iT+kk*6,iT+kk*6-6,j] <- FOI_eq[i,j]
        A[[i]][iP+kk*6,iP+kk*6-6,j] <- FOI_eq[i,j]
      }
      
    }
  }
  
  # Calculate equilibrium population in each compartment at age 0:
  for(j in 1:nh) {
    b[[1]]  <-  array(c(-p$eta*het_wt[j],rep(0,n-1)), dim=c(1,n,nh))
    # Set dS0/dt at age 0 to births, all other differential equations to 0 (for equilibrium values)
    
    # Solve matrix:
    Z[1,j,] <- solve(A[[1]][,,j], b[[1]][,,j])
  }
  
  # For all other age groups:
  for(j in 1:nh) {
    for(i in 2:na) {
      b[[i]] <-  array(c(-delta[i]*Z[i-1,j,1:n]), dim=c(1,n,nh))
      # Set all differential equations to the proportion of the population aging in
      
      # Solve matrix:
      Z[i,j,] <- solve(A[[i]][,,j], b[[i]][,,j])
    }
  }
  
  # Fill in equilibrium states
  T_eq <- Z[,,iT+kk_vec*6]
  P_eq <- Z[,,iP+kk_vec*6]
  S_eq <- Z[,,iS+kk_vec*6]
  I_PCR_eq <- Z[,,iI_PCR+kk_vec*6]
  I_LM_eq <- Z[,,iI_LM+kk_vec*6]
  I_D_eq <- Z[,,iI_D+kk_vec*6]
  
  
  ##############################
  ##############################
  # Vector model
  w_mean_blood_meal_rate <- sum(p$blood_meal_rates * p$species_proportions)
  w_mean_Q0 <- sum(p$Q0 * p$species_proportions)
  
  FOIvij_eq <- array(dim=c(na,nh,(K_max+1)))
  for (kk in 1:(K_max+1))
  {
    for (j in 1:nh)
    {
      for (i in 1:na)
      {
        FOIvij_eq[i, j, kk] <-  w_mean_blood_meal_rate * w_mean_Q0 * foi_age[i] *
          rel_foi[j]/omega * (p$ct * T_eq[i, j,kk] +
                                p$cd *I_D_eq[i, j,kk] + 
                                p$ca * I_LM_eq[i, j,kk] +
                                p$cu * I_PCR_eq[i, j,kk])

      }
    }
  }
  # Mosquito states
  FOIv_eq <- sum(FOIvij_eq)
  # Iv_eq <- FOIv_eq * p$Surv0/(FOIv_eq + p$mu0)
  # Sv_eq <- p$mu0 * Iv_eq/(FOIv_eq * p$Surv0)
  # Ev_eq <- 1 - Sv_eq - Iv_eq
  
  # Mosquito density needed to give this EIR
  # mv0 <- omega * EIRd_eq/(Iv_eq * p$av0)
  
  # Larval states
  # K0 <- 2 * mv0 * p$dLL * p$mu0 * (1 + p$dPL * p$muPL) * p$gammaL * (p$lambda + 1)/(p$lambda/(p$muLL *
  #                                                                                                             p$dEL) - 1/(p$muLL * p$dLL) - 1)
  # PL_eq <- 2 * p$dPL * p$mu0 * mv0
  # LL_eq <- p$dLL * (p$muPL + 1/p$dPL) * PL_eq
  # EL_eq <- (LL_eq/p$dLL + p$muLL* LL_eq * (1 + p$gammaL * LL_eq/K0))/(1/p$dEL - p$muLL * p$gammaL * LL_eq/K0)
  
  # Add in final dimension - interventions
  # num_int <- p$num_int
  # cov <- p$cov
  
  # mat <- array(0, dim=c(na,nh,nk))
  
  # S_eq <- vapply(cov, FUN = function(x)
  # {
  #   x * S_eq
  # }, mat)
  # T_eq <- vapply(cov, FUN = function(x)
  # {
  #   x * T_eq
  # }, mat)
  # I_D_eq <- vapply(cov, FUN = function(x)
  # {
  #   x * I_D_eq
  # }, mat)
  # I_LM_eq <- vapply(cov, FUN = function(x)
  # {
  #   x * I_LM_eq
  # }, mat)
  # I_PCR_eq <- vapply(cov, FUN = function(x)
  # {
  #   x * I_PCR_eq
  # }, mat)
  # P_eq <- vapply(cov, FUN = function(x)
  # {
  #   x * P_eq
  # }, mat)
  # Hyp_eq <- vapply(cov, FUN = function(x)
  # {
  #   x * Hyp_eq
  # }, mat)
  
  AP_eq = array(AP_eq, c(na, nh, nk))
  AC_eq = array(AC_eq, c(na, nh, nk))
  AP_MAT_eq = array(AP_MAT_eq, c(na, nh, nk))
  AC_MAT_eq = array(AC_MAT_eq, c(na, nh, nk))
  
  # Change array structure to: [age, het groups, intervention groups, hypnozoites]
  S_eq <- aperm(S_eq, c(1,2,3))
  I_PCR_eq <- aperm(I_PCR_eq, c(1,2,3))
  I_LM_eq <- aperm(I_LM_eq, c(1,2,3))
  I_D_eq <- aperm(I_D_eq, c(1,2,3))
  T_eq <- aperm(T_eq, c(1,2,3))
  P_eq <- aperm(P_eq, c(1,2,3))
  Hyp_eq <- aperm(Hyp_eq, c(1,2,3))
  AP_eq <- aperm(AP_eq, c(1,2,3))
  AC_eq <- aperm(AC_eq, c(1,2,3))
  AP_MAT_eq <- aperm(AP_MAT_eq, c(1,2,3))
  AC_MAT_eq <- aperm(AC_MAT_eq, c(1,2,3))
  
  
  # # Seasonality
  # admin_units_seasonal <- load_file("admin_units_seasonal.rds")
  # admin_matches <- admin_match(admin_unit = admin_unit, country = country,
  #                              admin_units_seasonal = admin_units_seasonal)
  
  # if(admin_matches == 0){
  #   ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
  # } else {
  #   ssa0 <- admin_units_seasonal$a0[admin_matches]
  #   ssa1 <- admin_units_seasonal$a1[admin_matches]
  #   ssa2 <- admin_units_seasonal$a2[admin_matches]
  #   ssa3 <- admin_units_seasonal$a3[admin_matches]
  #   ssb1 <- admin_units_seasonal$b1[admin_matches]
  #   ssb2 <- admin_units_seasonal$b2[admin_matches]
  #   ssb3 <- admin_units_seasonal$b3[admin_matches]
  #   theta_c <- admin_units_seasonal$theta_c[admin_matches]
  # }
  
  # better het bounds for equilbirum initialisation in individual model
  # zetas <- rlnorm(n = 1e5,meanlog = -p$sigma_squared/2, sdlog = sqrt(p$sigma_squared))
  # while(sum(zetas>100)>0){
  #   zetas[zetas>100] <- rlnorm(n = sum(zetas>100),meanlog = -p$sigma_squared/2, sdlog = sqrt(p$sigma_squared))
  # }
  
  # wt_cuts <- round(cumsum(het_wt)*1e5)
  # zeros <- which(wt_cuts==0)
  # wt_cuts[zeros] <- 1:length(zeros)
  # larges <- which(wt_cuts==1e5)
  # wt_cuts[larges] <- (1e5 - (length(larges)-1)):1e5
  # wt_cuts <- c(0,wt_cuts)
  # het_bounds <- sort(zetas)[wt_cuts]
  # het_bounds[length(het_bounds)] <- (max_age/365)+1
  
  # browser()
  
  c(sum(S_eq),
    sum(I_D_eq),
    sum(I_LM_eq),
    sum(I_PCR_eq),
    sum(T_eq),
    sum(P_eq))
  
  sum(0:10*colSums(colSums(Hyp_eq)))
  
  
  sum(rowSums(AP_eq*hyp_wt, dims = c(2))*outer(den,het_wt))
  sum(rowSums(AC_eq*hyp_wt, dims = c(2))*outer(den,het_wt))
  
  } else {
    
    # calculate immunity functions and onward infectiousness at equilibrium for
    # all age groups. See doi:10.1186/s12936-016-1437-9 for details of derivation.
    IC <- ID <- Hyp <- 0
    IDA <- ICA <- FOI <-FOIH <- q <- cA <- EPS <- HH <- rep(0, n_age)
    for (i in 1:n_age) {
      
      # rate of ageing plus death
      re <- r[i] + p$eta
      
      # 
      eps <- EIR*psi[i]
      EPS[i] <- eps
      
      # calculate probability of infection 
      b <- p$b
      
      # calculate force of infection (lambda)
      FOI[i] <- b*eps
      
      # calculate hypnozoite batch number
      Hyp <- (FOI[i] + re * Hyp) / (p$gammal + re)
      Hyp <- ifelse(Hyp>10, 10, Hyp)
      HH[i] <- Hyp
      
      # calculate force of infection with hypnozoites (lambda_H)
      FOIH[i] <- FOI[i] + HH[i]*p$f
      
      # update clinical immunity IC
      IC <- (FOIH[i]/(FOIH[i]*p$uc + 1) + re*IC)/(1/p$rc + re)
      ICA[i] <- IC
      
      # update detection immunity ID
      ID <- (FOIH[i]/(FOIH[i]*p$ud + 1) + re*ID)/(1/p$rid + re)
      IDA[i] <- ID
      
    }
    
    # calculate maternal clinical and severe immunity, assumed at birth to be a
    # proportion of the acquired immunity of a 20 year old
    ICM0 <- ICA[age20]*p$pcm
    IDM0 <- IDA[age20]*p$pcm
    
    ICM <- rep(0, n_age)
    IDM <- rep(0, n_age)
    for (i in 1:n_age) {
      # maternal clinical and severe immunity decays from birth
      if (i == n_age) {
        ICM[i] <- 0
        IDM[i] <- 0
      } else {
        ICM[i] <- ICM0 * p$rm / (age_days[i + 1] - age_days[i]) * (exp(-age_days[i] / (p$rm)) - exp(-age_days[i + 1] / (p$rm)))
        IDM[i] <- IDM0 * p$rm / (age_days[i + 1] - age_days[i]) * (exp(-age_days[i] / (p$rm)) - exp(-age_days[i + 1] / (p$rm)))
      }
    }
    
    # calculate probability of acquiring patent infection as a function of
    # different immunity types
    phi_patent <- p$philm_min + (p$philm_max-p$philm_min)/(1 + ((IDA+IDM)/p$alm50)^p$klm)
    
    # calculate probability of acquiring clinical disease as a function of
    # different immunity types
    phi_clin <- p$phi0*p$phi1  + (p$phi0 - p$phi0*p$phi1)/(1 + ((ICA+ICM)/p$ic0)^p$kc)
    
    # rate of subpatent recovery
    r_PCR  = 1/(p$dpcr_min + (p$dpcr_max - p$dpcr_min)/(1 + ((IDA+IDM)/p$apcr50)^p$kpcr))
    
    # calculate equilibrium solution of all model states. Again, see
    # doi:10.1186/s12936-016-1437-9 for details
    pos_LM <- pos_PCR <- inc <- pat_inc <- clin_inc <- rep(0, n_age)
    S <- T <- P <- D <- A <- U <- Y <- rep(0, n_age)
    
    
    MM_i <- function(i)
    {
      MM <- matrix(0, nrow=6, ncol=6)
      
      MM[1,1] = - FOIH[i] - p$eta - r[i]
      MM[1,2] = + r_PCR[i]
      MM[1,6] = + p$rp
      
      MM[2,1] = + FOIH[i]*(1-phi_patent[i])
      MM[2,2] = + FOIH[i]*(1-phi_patent[i]) - FOIH[i] - r_PCR[i] - p$eta - r[i]
      MM[2,3] = + 1/p$da
      
      MM[3,1] = + FOIH[i]*phi_patent[i]*(1-phi_clin[i])
      MM[3,2] = + FOIH[i]*phi_patent[i]*(1-phi_clin[i])
      MM[3,3] = + FOIH[i]*(1-phi_clin[i]) - FOIH[i] - 1/p$da - p$eta - r[i]
      MM[3,4] = + 1/p$dd
      
      MM[4,1] = + FOIH[i]*phi_patent[i]*phi_clin[i]*(1-ft)
      MM[4,2] = + FOIH[i]*phi_patent[i]*phi_clin[i]*(1-ft)
      MM[4,3] = + FOIH[i]*phi_clin[i]*(1-ft)
      MM[4,4] = - 1/p$dd - p$eta - r[i]
      
      MM[5,1] = + FOIH[i]*phi_patent[i]*phi_clin[i]*ft
      MM[5,2] = + FOIH[i]*phi_patent[i]*phi_clin[i]*ft
      MM[5,3] = + FOIH[i]*phi_clin[i]*ft
      MM[5,5] = - 1/p$dt - p$eta - r[i]
      
      MM[6,5] = + 1/p$dt
      MM[6,6] = - p$rp - p$eta - r[i]
      
      MM
    }
    
    S <- U <- A <-D <- T <- P <- rep(NA, n_age)
    
    MM <- MM_i(1)
    BB <- rep(0, 6)
    BB[1] = - p$eta
    
    XX <- solve(MM)%*%BB
    
    S[1] <- XX[1]
    U[1] <- XX[2]
    A[1] <- XX[3]
    D[1] <- XX[4]
    T[1] <- XX[5]
    P[1] <- XX[6]
    
    for(i in 2:n_age)
    {
      MM <- MM_i(i)
      BB = - r[i-1]*XX
      
      XX = solve(MM)%*%BB
      
      S[i] <- XX[1]
      U[i] <- XX[2]
      A[i] <- XX[3]
      D[i] <- XX[4]
      T[i] <- XX[5]
      P[i] <- XX[6]
      
      pos_LM[i] <- D[i] + T[i] + A[i]
      pos_PCR[i] <- D[i] + T[i] + A[i] + U[i]
      
      # calculate clinical incidence
      Y[i] <- S[i]+U[i]+A[i]
      inc[i] <- Y[i]*FOIH[i]
      pat_inc[i] <- Y[i]*FOIH[i]*phi_patent[i]
      clin_inc[i] <- Y[i]*FOIH[i]*phi_patent[i]*phi_clin[i]
    }
    
    # calculate mean infectivity
    inf <- p$cd*D + p$ct*T + p$ca*A + p$cu*U
    
    
    
    
  }
  
  ## collate init
  state <- list(init_S = S_eq, init_T = T_eq, init_I_D = I_D_eq, init_I_LM = I_LM_eq,
                init_I_PCR = I_PCR_eq, init_P = P_eq, init_Hyp = Hyp_eq,
                init_AP = AP_eq, init_AC = AC_eq, init_AP_MAT = AP_MAT_eq,
                init_AC_MAT = AC_MAT_eq, #init_Iv = Iv_eq, init_Sv = Sv_eq,
                # init_Ev = Ev_eq, init_PL = PL_eq, init_LL = LL_eq, init_EL = EL_eq,
                age_width = age_width, age_rate = age_rate, het_wt = het_wt, het_x = het_x,
                omega = omega, foi_age = foi_age, rel_foi = rel_foi,
                #K0 = K0, mv0 = mv0, 
                na = na, nh = nh,
                nk=nk, kk=kk_vec, K_max = K_max, K_max_switch = K_max_switch,
                K0_switch = K0_switch, K_max_switch_on = K_max_switch_on,
                #ni = num_int,
                hyp_wt = hyp_wt,
                FOI = FOI_eq, EIR_eq = EIR_eq, Hyp_eq = Hyp_eq,
                phi_LM_eq = phi_LM_eq, phi_D_eq = phi_D_eq, dPCR_eq = dPCR_eq,
                den = den, #age59 = age59, age05 = age05, ssa0 = ssa0, ssa1 = ssa1,
                # ssa2 = ssa2, ssa3 = ssa3, ssb1 = ssb1, ssb2 = ssb2, ssb3 = ssb3,
                # theta_c = theta_c, 
                age = age, ft = ft, FOIv_eq = FOIv_eq,
                # FOIvij_eq=FOIvij_eq,
                age_mid_point = age_mid_point, #het_bounds = het_bounds,
                pi = pi,
                age20l = age20l, age20u = age20u, age_20_factor = age_20_factor)
  
  state <- append(state,p)
  
  # browser()
  ## Prep results for malariasimultion:
  if(isTRUE(malariasimulationoutput)){
    ## List based on het group
    ## Summarise over hypnozoites
    malsimres <- lapply(1:nh, function(h){
      data.frame("age" = age[-length(age)]/365,
                 "S" = rowSums(S_eq[,h,]),
                 "D" = rowSums(I_D_eq[,h,]),
                 "A" = rowSums(I_LM_eq[,h,]),
                 "U" = rowSums(I_PCR_eq[,h,]),
                 "T" = rowSums(T_eq[,h,]),
                 "P" = rowSums(P_eq[,h,]),
                 "ICA" = rowSums(AC_eq[,h,] * hyp_wt[,h,]),
                 "ID" = rowSums(AP_eq[,h,] * hyp_wt[,h,]),
                 "ICM" = rowSums(AC_MAT_eq[,h,] * hyp_wt[,h,]),
                 "IDM" = rowSums(AP_MAT_eq[,h,] * hyp_wt[,h,]),
                 "HH" = rowSums(t(t(hyp_wt[,h,])*0:K_max)),
                 "psi" = foi_age,
                 "prop" = den,
                 "phi_clin" = rowSums(phi_D_eq[,h,] * hyp_wt[,h,]))
                 # "inf" = rowSums(inf[,h,] * hyp_wt[,h,]))
    })
    return(list(ret = malsimres, FOIM = FOIv_eq))
  }
  
  return(state)
  
}
