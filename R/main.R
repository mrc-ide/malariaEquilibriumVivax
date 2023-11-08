#------------------------------------------------
#' @title Vivax equilibrium solution with biting heterogeneity
#'
#' @description Returns the vivax equilibrium states for the model of White et al.,
#'   (2018).
#'
#'   This function does account for biting heterogeneity and return the 
#'   equilibria listed over heterogeneity - to summarise and calculate FOIM 
#'   please use \code{human_equilibrium_vivax_summarise_het()}.
#'   
#'   This version uses Michael White's code as a template. I've removed the omega 
#'   standardisation for consistency with other equiilbria.
#'
#' @param EIR EIR for adults, in units of infectious bites per person per year
#' @param ft proportion of clinical cases effectively treated
#' @param p vector of model parameters
#' @param age vector of age groups, in units of years
#' @param h a list of Gauss-Hermite nodes and associated weights, used for
#'   integrating over heterogeneity in biting. See \code{?gq_normal} for an
#'   example.
#'
#' @export
#' 
human_equilibrium_vivax_full_het <- function(EIR, ft, p, age, h = NULL, no_omega = NULL, use_simple_ages = NULL, K_max = 10) {
  
  # This code is based on Michael White code
  # tailored to malariasimulation parameterisation.
  
  # check inputs
  assert_single_pos(EIR, zero_allowed = FALSE)
  assert_single_bounded(ft)
  assert_custom_class(p, "list")
  assert_vector_pos(age)
  assert_noduplicates(age)
  assert_increasing(age)
  
  ###################################
  ## 1. ##  Age and heterogeneity  ##
  ###################################
  
  #####################################
  ## Age demography
  # some changes required to MW code when using the EQUILIBRIUMAGES or
  # equivalent input (age) from malariasimulation
  
  # get basic properties
  n_age <- length(age)
  age_days <- age*365
  EIR <- EIR/365
  
  # parameter conversions
  eta <- 1/p$average_age
  r_par <- 1/p$rid
  r_clin <- 1/p$rc
  rp <- 0.1
  rd <- 1/p$dd
  rt <- 1/p$dt
  ra <- 1/p$da
  ru <- 1/p$du
  
  
  if(is.null(use_simple_ages)){
    ## Calculate prop: the proportion within each age bracket
    ## These methods are slightly different from malariaEquilibrium
    ## but produce similar estimates
    ## prop is equivalent to prop in MW equilibrium
    ## age_days are eqiovalent to age_bounds, without the upper bound
    ## mean_age is 1/p$eta
    prop <- rep(NA, n_age)
    prop[1:(n_age-1)]  <-  exp( - age_days[1:(n_age-1)]*eta) - exp( - age_days[2:n_age]*eta )
    prop[n_age]        <- 1 - sum( prop[1:(n_age-1)] )
    
    ## Calculate rate of aging, 1/duration in age group
    ## r is equivalent to r in MW eq
    r <- rep(NA, n_age)
    for(i in 1:(n_age-1)){
      r[i] <- ( 1 - sum(prop[1:i]) )/( prop[i]/eta )
    }
    
    r[n_age] = 0
    
    ###################################
    ##
    ## Modification of aging rates for immune functions.
    ##
    ## Proportions infected must always accord to the demographic
    ## constraints when stratified. Immunity is a scalar and can
    ##(in theory) increase without bound. We must make an adjustment
    ## due to demography for aging of immunity.
    ## Provide better explanation.
    
    r_imm_age <- rep(NA, n_age)
    r_imm_age[1:(n_age-1)] <- r[1:(n_age-1)]*(prop[1:(n_age-1)]/prop[2:n_age])
    r_imm_age[n_age] <- 0
    
  } else {
    prop <- r <- rep(0, n_age)
    for (i in 1:n_age) {
      
      # r[i] can be thought of as the rate of ageing in this age group, i.e.
      # 1/r[i] is the duration of this group
      if (i == n_age) {
        r[i] <- 0
      } else {
        age_width <- age_days[i+1] - age_days[i]
        r[i] <- 1/age_width
      }
      
      # prop is calculated from the relative flows into vs. out of this age group.
      # For the first age group the flow in is equal to the birth rate (eta). For
      # all subsequent age groups the flow in represents ageing from the previous
      # group. The flow out is always equal to the rate of ageing or death.
      if (i == 1) {
        prop[i] <- eta/(r[i] + eta)
      } else {
        prop[i] <- prop[i-1]*r[i-1]/(r[i] + eta)
      }
    }
    r_imm_age <- r + eta
  }
  
  # calculate midpoint of age range. There is no midpoint for the final age
  # group, so use beginning of range instead
  ## age_days_midpoint is age_mids in MW
  age_days_midpoint = c(0.5*( age_days[2:(n_age)] + age_days[1:(n_age-1)]), age_days[n_age])
  index_MI_20 <- which.min((age_days_midpoint - 20*365)^2)[1]   ## index for the age group of a 20 year old (for maternal immunity)
  
  # calculate relative biting rate for each age group
  ## p$rho is rho_age
  ## p$a0 os age_0
  ## psi is age_bite
  psi <- 1 - p$rho*exp(-age_days_midpoint/p$a0)
  # Here we factor psi (age-specific biting rate) by omega (age-normalisation)
  # In malariaEquilibrium this stage is done later in the human_equilibrium function
  # not in the human_equilibrium_no_het function
  omega_age <- 1/sum(prop * psi)
  
  ### Nora's code:
  if(is.null(no_omega)){
    psi <- omega_age*psi  
  }
  ###
  
  ###########################################################################
  ## Heterogeneity in mosquito bites
  # In malariaEquilibrium this stage is done later in the human_equilibrium function
  # not in the human_equilibrium_no_het function
  if(isTRUE(p$enable_heterogeneity)){
    sig_het = sqrt(p$sigma_squared)
    x_het <- exp(-p$sigma_squared*0.5 + sqrt(p$sigma_squared)*h$nodes)
    w_het <- h$weights
    x_age_het <- psi%o%x_het
    w_age_het <- prop%o%w_het
    n_het <- length(h$nodes)
  } else {
    x_age_het <- psi
    w_het <- 1
    n_het <- 1
  }
  
  
  #################################################
  ## 2. ##  FoI and proportion with hypnozoites  ##
  #################################################
  
  lam_eq = as.matrix(EIR*p$b*x_age_het)
  ###### Nora's code:
  # lam_eq = as.matrix(EIR*p$b*x_age_het)/sum(x_het*w_het)
  ######
  HH_eq <- matrix(NA, nrow=n_age, ncol=n_het)
  for(j in 1:n_het){
    HH_eq[1,j] = ( 1/(p$gammal + r_imm_age[1]) )*( lam_eq[1,j] )
    
    for(i in 2:n_age){
      HH_eq[i,j] = ( 1/(p$gammal + r_imm_age[i]) )*( lam_eq[i,j] + r_imm_age[i]*HH_eq[i-1,j] )
    }
  }
  
  #######################
  # Should we cap the batch number?
  HH_eq[HH_eq>K_max] <- K_max
  #######################
  
  # Is this equivalent to modelling infections as two processes?
  lam_H_eq = lam_eq + p$f*HH_eq
  
  #################################################
  ## 3. ##  Equilibrium levels of immunity       ##
  #################################################
  
  A_par_eq <- matrix(NA, nrow=n_age, ncol=n_het)
  A_clin_eq <- matrix(NA, nrow=n_age, ncol=n_het)
  
  for(j in 1:n_het){
    A_par_eq[1,j] = ( 1/(r_par + r_imm_age[1]) )*( lam_H_eq[1,j]/(lam_H_eq[1,j]*p$ud+1) )
    
    for(i in 2:n_age){
      A_par_eq[i,j] = ( 1/(r_par + r_imm_age[i]) )*( lam_H_eq[i,j]/(lam_H_eq[i,j]*p$ud+1) + r_imm_age[i]*A_par_eq[i-1,j] )
    }
  }
  
  for(j in 1:n_het){
    A_clin_eq[1,j] = ( 1/(r_clin + r_imm_age[1]) )*( lam_H_eq[1,j]/(lam_H_eq[1,j]*p$uc+1) )
    
    for(i in 2:n_age){
      A_clin_eq[i,j] = ( 1/(r_clin + r_imm_age[i]) )*( lam_H_eq[i,j]/(lam_H_eq[i,j]*p$uc+1) + r_imm_age[i]*A_clin_eq[i-1,j] )
    }
  }
  
  A_par_MI_eq  = (p$pcm*exp(-age_days_midpoint/p$rm) )%o%A_par_eq[index_MI_20,]
  A_clin_MI_eq = (p$pcm*exp(-age_days_midpoint/p$rm) )%o%A_clin_eq[index_MI_20,]
  
  #######################################
  ## Effects of immune functions
  
  r_PCR_eq  = 1/(p$dpcr_min + (p$dpcr_max - p$dpcr_min)/(1 + ((A_par_eq + A_par_MI_eq)/p$apcr50)^p$kpcr))
  phi_LM_eq = p$philm_min + (p$philm_max-p$philm_min)/(1 + ((A_par_eq + A_par_MI_eq)/p$alm50)^p$klm)
  phi_D_eq  = p$phi0*p$phi1  + (p$phi0 - p$phi0*p$phi1)/(1 + ((A_clin_eq + A_clin_MI_eq)/p$ic0)^p$kc)
  
  
  ##################################################
  ## 4. ##  Function for equilibrium solution     ##
  ##    ##  via linear algebra                    ##
  ##################################################
  
  MM_ij <- function(i, j)
  {
    MM <- matrix(0, nrow=6, ncol=6)

    MM[1,1] = - lam_H_eq[i,j] - eta - r[i]
    MM[1,2] = + r_PCR_eq[i,j]
    MM[1,6] = + rp
    
    MM[2,1] = + lam_H_eq[i,j]*(1.0-phi_LM_eq[i,j])
    MM[2,2] = + lam_H_eq[i,j]*(1.0-phi_LM_eq[i,j]) - lam_H_eq[i,j] - r_PCR_eq[i,j] - eta - r[i]
    MM[2,3] = + ra
    
    MM[3,1] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*(1.0-phi_D_eq[i,j])
    MM[3,2] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*(1.0-phi_D_eq[i,j])
    MM[3,3] = + lam_H_eq[i,j]*(1.0-phi_D_eq[i,j]) - lam_H_eq[i,j] - ra - eta - r[i]
    MM[3,4] = + rd
    
    MM[4,1] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*phi_D_eq[i,j]*(1-ft)
    MM[4,2] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*phi_D_eq[i,j]*(1-ft)
    MM[4,3] = + lam_H_eq[i,j]*phi_D_eq[i,j]*(1-ft)
    MM[4,4] = - rd - eta - r[i]
    
    MM[5,1] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*phi_D_eq[i,j]*ft
    MM[5,2] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*phi_D_eq[i,j]*ft
    MM[5,3] = + lam_H_eq[i,j]*phi_D_eq[i,j]*ft
    MM[5,5] = - rt - eta - r[i]
    
    MM[6,5] = + rt
    MM[6,6] = - rp - eta - r[i]
    
    MM
  }
  
  
  #################################################
  ## 5. ##  Solve for equilibrium solution       ##
  #################################################
  
  ###################################################
  ## Objects for storing equilibrium
  
  S_eq     <- matrix(NA, nrow=n_age, ncol=n_het)
  I_PCR_eq <- matrix(NA, nrow=n_age, ncol=n_het)
  I_LM_eq  <- matrix(NA, nrow=n_age, ncol=n_het)
  D_eq     <- matrix(NA, nrow=n_age, ncol=n_het)
  T_eq     <- matrix(NA, nrow=n_age, ncol=n_het)
  P_eq     <- matrix(NA, nrow=n_age, ncol=n_het)
  
  for(j in 1:n_het){
    
    
    MM <- MM_ij(1,j)
    
    BB    = rep(0, 6)
    BB[1] = - w_het[j]*eta
    
    XX = solve(MM)%*%BB
    
    S_eq[1,j]     = XX[1]
    I_PCR_eq[1,j] = XX[2]
    I_LM_eq[1,j]  = XX[3]
    D_eq[1,j]     = XX[4]
    T_eq[1,j]     = XX[5]
    P_eq[1,j]     = XX[6]
    
    
    for(i in 2:n_age)
    {
      MM <- MM_ij( i, j )
      
      BB = - r[i-1]*XX
      
      XX = solve(MM)%*%BB
      
      S_eq[i,j]     = XX[1]
      I_PCR_eq[i,j] = XX[2]
      I_LM_eq[i,j]  = XX[3]
      D_eq[i,j]     = XX[4]
      T_eq[i,j]     = XX[5]
      P_eq[i,j]     = XX[6]
    }
  }

  # Mean infectivity
  inf <- p$cd*D_eq + p$ct*T_eq + p$ca*I_LM_eq + p$cu*I_PCR_eq
  
  # calculate proportion detectable by mocroscopy and PCR
  # These previous had more complex equations to calculate A_eq for example, involving q and an aA parameter...?
  pos_M <- D_eq + T_eq + I_LM_eq
  pos_PCR <- D_eq + T_eq + I_LM_eq + I_PCR_eq
  
  # incidence
  inc <- (S_eq + I_LM_eq + I_PCR_eq) * lam_H_eq
  asym_inc <- inc * phi_LM_eq
  clin_inc <- asym_inc * phi_D_eq
  
  # ??? I'm not sure what EPS stands for...? (If someone could tell me, please let me know!)
  EPS <- psi*EIR
  
  
  ######################################
  ## Add back in omega to calculate FOIM
  # psi <- omega_age*psi
  ######################################
  
  # return matrix: Consider outputs?
  ret <- lapply(X = 1:n_het, FUN = function(j){
    cbind(
      age = age,
      HH = HH_eq[,j],
      S = S_eq[,j],
      T = T_eq[,j],
      D = D_eq[,j],
      A = I_LM_eq[,j],
      U = I_PCR_eq[,j],
      P = P_eq[,j],
      inf = inf[,j],
      prop = prop,
      psi = psi,
      pos_M = pos_M[,j],
      pos_PCR = pos_PCR[,j],
      inc = inc[,j],
      asym_inc = asym_inc[,j],
      clin_inc = clin_inc[,j],
      ICA = A_clin_eq[,j],
      ICM = A_clin_MI_eq[,j],
      ID = A_par_eq[,j],
      IDM = A_par_MI_eq[,j],
      FOI = lam_eq[,j],
      FOIH = lam_H_eq[,j],
      phi = phi_D_eq[,j],
      phi_lm = phi_LM_eq[,j],
      rpcr = r_PCR_eq[,j],
      EPS = EPS,
      r = r,
      r_imm_age = r_imm_age
    )
  })

  if(isTRUE(p$enable_heterogeneity)){
    return(list(ret = ret, x_het = x_het, w_het = w_het))
    
    } else {
      return(list(ret = ret))
    }
  
}


#------------------------------------------------
#' @title Vivax equilibrium solution
#'
#' @description Returns the vivax equilibrium states for the model of White et al.,
#'   (2018).
#'
#' @param eq the equilibrium solution calcualted by human_equilibrium_vivax_full_het
#' @param p a parameter set generated by get_parameters() and complementary functions
#' in malariasimulation
#'
#' @export

human_equilibrium_vivax_summarise <- function(eq, p){
  
  #######################################################
  ## 1. ##  Age stratify                               ##
  #######################################################

  ## sum states
  states <- Reduce("+", lapply(eq$ret, function(x){x[,c("S","T","D","A","U","P","inf")]}))
  ## age
  age <- eq$ret[[1]][,1]
  
  if(isTRUE(p$enable_heterogeneity)){
    x_het <- eq$x_het
    w_het <- eq$w_het
    
    ## Weight hypnozoites and immunities
    hh_imm <- Reduce("+",
                     lapply(1:length(w_het),function(x){eq$ret[[x]][,c("HH","ICA","ICM","ID","IDM","phi","prop","FOI","FOIH")] * w_het[x]}))
    # foim <- sum(unlist(lapply(1:length(w_het),function(x){sum(eq$ret[[x]][,c("inf")] * eq$ret[[x]][,c("psi")]) * w_het[x] * x_het[x]})))
    foim <- sum(unlist(lapply(1:length(w_het),function(x){sum(eq$ret[[x]][,c("inf")] * eq$ret[[x]][,c("psi")]) * x_het[x]})))
  } else {
    hh_imm <- eq$ret[[1]][,c("HH","ICA","ICM","ID","IDM","phi","prop")]
    foim <- sum(eq$ret[[1]][,c("inf")] * eq$ret[[1]][,c("psi")])
  }
  
  #######################################################
  ## 2. ##  FOIM                                       ##
  #######################################################
  
  # FOIM is foim in MW equilibrium
  # zeta (heteogeneity) and omega (normalising over age) are taken care of within the equilibrium solution
  # the only thing that remains is to multiply by alpha: rate at which a mosquito takes a blood meal on humans
  
  ## Get ages
  # a <- parameters$
  
  ## get psi
  # psi <- states[,"psi"]
  ## get zeta
  
  ## get prob bitten
  # prob_bitten <- 1
  ## get .pi?
  ## get Q0
  # p$Q0
  ## get W
  ## get Z
  ## Get f
  # p$blood_meal_rates
  
  # hh_imm[,"psi"]
  
  # weighted_infectivity <- sum(colSums(states[,c("T","D","A","U")]) * hh_imm[,"psi"] * c(unlist(p[c("ct","cd","ca","cu")])))
  eta <- 1/p$average_age
  if(!is.null(no_omega)){
    omega <- 1 - p$rho*eta/(eta + 1/p$a0)
  } else {
    omega <- 1
  }
  
  alpha <- p$blood_meal_rates * p$Q0
  foim <- foim * alpha / omega
  eq_summary <- cbind(age, states, hh_imm)
  
  return(list(states = eq_summary, FOIM = foim))
}
