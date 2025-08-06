#------------------------------------------------
#' @title Vivax equilibrium solution with biting heterogeneity
#'
#' @description Returns the vivax equilibrium states for the model of White et al.,
#'   (2018).
#'
#' @param EIR EIR for adults, in units of infectious bites per person per year
#' @param ft proportion of clinical cases effectively treated
#' @param p vector of model parameters
#' @param age vector of age groups, in units of years
#'
#' @export

vivax_equilibrium <- function(EIR, ft, p, 
                              age = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
                                      1.2, 1.4, 1.6, 1.8, 2.0,
                                      2.2, 2.4, 2.6, 2.8, 3.0,
                                      3.4, 3.8, 4.2, 4.6, 5.0,
                                      5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
                                      11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                                      22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0,
                                      45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0)){
  
  ## Check Parameters
  if(!is.numeric(age)) stop("age provided is not numeric")
  if(!is.numeric(p$n_heterogeneity_groups)) stop("het_brackets provided is not numeric")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  if(!is.numeric(EIR)) stop("EIR provided is not numeric")
  if(!is.numeric(p$kmax)) stop("K_max is not provided or not numeric")
  if(!(p$kmax %in% c(1:30))) stop("K_max is not an integer between 1 and 30")
  
  assert_single_pos(EIR, zero_allowed = FALSE)
  assert_single_bounded(ft)
  assert_custom_class(p, "list")
  assert_vector_pos(age)
  assert_noduplicates(age)
  assert_increasing(age)
  
  #####################################
  ## 1.1. ##  Age and heterogeneity  ##
  #####################################  
  
  ##################################### 
  ## Age demography - replace with our own age demography, acknowledging that the oldest age group is dropped
  ## These represent age bounds
  max_age  = max(age)*365
  mean_age <- p$mean_age
  age_0 <- p$age_0
  sigma_squared <- p$sigma_squared
  rho_age <- p$rho_age
  
  mu_H <- 1/mean_age
  age_bounds <- 365 * age
  N_age <- length(age_bounds) - 1
  
  age_width = age_bounds[2:(N_age+1)] - age_bounds[1:N_age]
  age_mids = 0.5*( age_bounds[2:(N_age+1)] + age_bounds[1:N_age] )
  index_MI_20 <- which.min( (age_mids - 20*365)^2 )[1]   ## index for the age group of a 20 year old (for maternal immunity)
  
  age_demog <- rep(NA, N_age)
  age_demog[1:(N_age-1)]  <-  exp( - age_bounds[1:(N_age-1)]/mean_age ) - exp( - age_bounds[2:N_age]/mean_age ) 
  age_demog[N_age]        <- 1 - sum( age_demog[1:(N_age-1)] )
  
  r_age <- rep(NA, N_age)
  for(i in 1:(N_age-1)){
    r_age[i] <- ( 1 - sum(age_demog[1:i]) )/( age_demog[i]*mean_age )
  }	
  r_age[N_age] = 0
  
  ###########################################################################
  ## Age-dependent exposure to mosquito bites
  
  age_bite = 1 - rho_age*exp(-age_mids/age_0)
  omega_age = 1/sum(age_demog*age_bite)
  ## This part standardises the population level (including age structuring) to
  ## the initial EIR value. The Pf model does not do this, so for consistency we
  ## are changing this to the Pf method. We do not believe this will make a 
  ## fundemental change to analytical outputs.
  if(isTRUE(p$vivax_EIR_at_population_level)){
    age_bite = omega_age*age_bite
  }
  
  ###########################################################################
  ## Heterogeneity in mosquito bites
  
  sig_het = sqrt(sigma_squared)
  N_het <- p$n_heterogeneity_groups
  x_het <- exp(statmod::gauss.quad.prob(N_het, dist="normal", mu=-0.5*sig_het^2, sigma=sig_het)$nodes)
  w_het <-     statmod::gauss.quad.prob(N_het, dist="normal", mu=-0.5*sig_het^2, sigma=sig_het)$weights
  
  x_age_het <- age_bite%o%x_het
  w_age_het <- age_demog%o%w_het
  
  ################################################################
  ## 1.2. ##  Parameter definitions and additional assignments  ##
  ################################################################
  
  aa            <- p$Q0 * p$blood_meal_rates
  bb            <- p$bb           ## mosquito -> human transmission probability
  
  c_PCR         <- p$c_PCR        ## human -> mosquito transmission probability (PCR)
  c_LM          <- p$c_LM         ## human -> mosquito transmission probability (LM-detectable)
  c_D           <- p$c_D          ## human -> mosquito transmission probability (disease state)
  c_T           <- p$c_T          ## human -> mosquito transmission probability (treatment)
  
  d_E           <- p$d_E          ## duration of liver-stage latency
  r_D           <- p$r_D          ## duraton of disease = 1/rate
  r_T           <- p$r_T          ## duraton of prophylaxis = 1/rate
  r_P           <- 1/10           ## duraton of treatment = 1/rate
  chi_treat     <- ft             ## proportion of symptomatic episodes receiving first-line treatment
  
  r_par         <- p$r_par        ## rate of decay of anti-parasite immunity
  r_clin        <- p$r_clin       ## rate of decay of clinical immunity
  
  
  d_PCR_min     <- p$d_PCR_min    ## maximum duration of PCR-detectable infection
  
  mu_M          <- p$mu_M         ## mosquito death rate = 1/(mosquito life expectancy): 1/6
  tau_M         <- p$tau_M        ## duration of sporogony
  
  ff            <- p$ff           ## relapse rate: 1/41
  gamma_L       <- p$gamma_L      ## duration of liver-stage carriage: 1/383
  
  K_max         <- p$K_max        ## Maximum number of hypnozoite batches
  
  u_par         <- p$u_par        ## refractory period for anti-parasite immune boosting
  phi_LM_max    <- p$phi_LM_max   ## probability of LM_detectable infection with no immunity
  phi_LM_min    <- p$phi_LM_min   ## probability of LM_detectable infection with maximum immunity
  A_LM_50pc     <- p$A_LM_50pc    ## blood-stage immunity scale parameter
  K_LM          <- p$K_LM         ## blood-stage immunity shape parameter
  u_clin        <- p$u_clin       ## refractory period for clinical immune boosting
  phi_D_max     <- p$phi_D_max    ## probability of clinical episode with no immunity
  phi_D_min     <- p$phi_D_min    ## probability of clinical episode with maximum immunity
  A_D_50pc      <- p$A_D_50pc     ## clinical immunity scale parameter
  K_D           <- p$K_D          ## clinical immunity shape parameter
  A_d_PCR_50pc  <- p$A_d_PCR_50pc ## scale parameter for effect of anti-parasite immunity on PCR-detectable infection
  K_d_PCR       <- p$K_d_PCR      ## shape parameter for effect of anti-parasite immunity on PCR-detectable infection
  d_PCR_max     <- p$d_PCR_max    ## maximum duration on PCR-detectable infection
  d_LM          <- p$d_LM         ## duration of LM-detectable infection
  P_MI          <- p$P_MI         ## Proportion of immunity acquired maternally
  d_MI          <- p$d_MI         ## Rate of waning of maternal immunity
  
  r_LM <- 1/d_LM
  
  EIR_site   = EIR/365
  
  
  #######################################
  ## 1.3. ##  Equilibrium calculation  ##
  #######################################
  
  K_MAT = matrix(0, nrow=K_max+1, ncol=K_max+1)
  diag(K_MAT) = 0:K_max
  
  D_MAT = diag(K_max+1)
  
  L_MAT = matrix(0, nrow=K_max+1, ncol=K_max+1)
  for(k in 1:K_max)
  {
    L_MAT[k,  k+1] = + k
    L_MAT[k+1,k+1] = - k
  }
  
  OD_MAT = matrix(0, nrow=K_max+1, ncol=K_max+1)
  for(k in 1:K_max)
  {
    OD_MAT[k+1,k] = + 1
  }
  OD_MAT[K_max+1,K_max+1] = 1
  
  
  H_MAT = matrix(0, nrow=K_max+1, ncol=K_max+1)
  for(k in 1:K_max)
  {
    H_MAT[k,k]   = - 1
    H_MAT[k+1,k] = + 1
  }
  
  S_index     = 0*(K_max+1)+(1:(K_max+1))
  I_PCR_index = 1*(K_max+1)+(1:(K_max+1))
  I_LM_index  = 2*(K_max+1)+(1:(K_max+1))
  I_D_index   = 3*(K_max+1)+(1:(K_max+1))
  T_index     = 4*(K_max+1)+(1:(K_max+1))
  P_index     = 5*(K_max+1)+(1:(K_max+1))
  
  
  ###################################################
  ## 2.1. ##  FoI and proportion with hypnozoites  ##
  ###################################################
  
  lam_eq = EIR_site*bb*x_age_het
  
  HH_eq = array(NA, dim=c(N_age, N_het, K_max+1))
  for(j in 1:N_het){
    ##########################
    ## Youngest age category
    
    HH_vec = rep(0, K_max+1)
    HH_vec[1] = w_het[j]
    HH_ik = lam_eq[1,j]*( H_MAT ) + gamma_L*( L_MAT ) - mu_H*( D_MAT ) - r_age[1]*( D_MAT ) 
    HH_vec = - solve(HH_ik)%*%(mu_H*HH_vec)
    HH_eq[1,j,] = HH_vec
    
    ##########################
    ## Older age categories
    
    for(i in 2:N_age){
      HH_ik = lam_eq[i,j]*( H_MAT ) + gamma_L*( L_MAT ) - mu_H*( D_MAT ) - r_age[i]*( D_MAT ) 
      HH_vec = - solve(HH_ik)%*%(r_age[i-1]*HH_vec)
      HH_eq[i,j,] = HH_vec		
    }
  }
  
  
  lam_H_eq  = array(NA, dim=c(N_age, N_het, K_max+1))
  for(i in 1:N_age){
    for(j in 1:N_het){
      lam_H_eq[i,j,] = lam_eq[i,j] + ff*(0:K_max)
    }
  }
  
  
  ###########################################################
  ## 2.2. ##  Equilibrium levels of anti-parasite immunity ##
  ###########################################################
  
  A_par_eq  = array(NA, dim=c(N_age, N_het, K_max+1))
  A_par_mat_eq = array(NA, dim=c(N_age, N_het, K_max+1))
  A_par_total_eq = array(NA, dim=c(N_age, N_het, K_max+1))
  
  for(j in 1:N_het){
    #############################	
    ## Youngest age category
    
    G_VEC          = c(0, lam_eq[1,j]*HH_eq[1,j,1:K_max]/HH_eq[1,j,2:(K_max+1)] ) + (0:K_max)*ff
    G_VEC[K_max+1] = G_VEC[K_max+1] + lam_eq[1,j] 
    G_VEC          = G_VEC/( G_VEC*u_par + 1 )
    
    LAM_MAT = - diag(K_max+1)
    LAM_MAT[K_max+1,K_max+1] = 0
    for(k in 1:K_max){
      LAM_MAT[k+1,k] = HH_eq[1,j,k]/HH_eq[1,j,k+1]
    }
    
    GAM_MAT = - diag(0:K_max)
    for(k in 1:K_max){
      GAM_MAT[k,k+1] = k*HH_eq[1,j,k+1]/HH_eq[1,j,k]
    }
    
    Z_MAT = - r_par*D_MAT + lam_eq[1,j]*LAM_MAT + gamma_L*GAM_MAT - mu_H*D_MAT - r_age[1]*D_MAT
    A_par_vec = solve( -Z_MAT )%*%( G_VEC )
    A_par_eq[1,j,] = A_par_vec
    
    ##########################
    ## Older age categories
    
    for(i in 2:N_age){
      G_VEC          = c(0, lam_eq[i,j]*HH_eq[i,j,1:K_max]/HH_eq[i,j,2:(K_max+1)] ) + (0:K_max)*ff
      G_VEC[K_max+1] = G_VEC[K_max+1] + lam_eq[i,j] 
      G_VEC          = G_VEC/( G_VEC*u_par + 1 )
      
      r_VEC = r_age[i-1]*A_par_eq[i-1,j,]*HH_eq[i-1,j,]/HH_eq[i,j,]
      
      LAM_MAT = - diag(K_max+1)
      LAM_MAT[K_max+1,K_max+1] = 0
      for(k in 1:K_max){
        LAM_MAT[k+1,k] = HH_eq[i,j,k]/HH_eq[i,j,k+1]
      }
      
      GAM_MAT = - diag(0:K_max)
      for(k in 1:K_max){
        GAM_MAT[k,k+1] = k*HH_eq[i,j,k+1]/HH_eq[i,j,k]
      }
      
      Z_MAT = - r_par*D_MAT + lam_eq[i,j]*LAM_MAT + gamma_L*GAM_MAT - mu_H*D_MAT - r_age[i]*D_MAT
      A_par_vec = solve( -Z_MAT )%*%( G_VEC + r_VEC )
      A_par_eq[i,j,] = as.vector(A_par_vec)
    }
  }
  
  
  ##########################
  ## Average over hypnozoites
  A_par_eq_ave = matrix(NA, nrow=N_age, ncol=N_het)
  
  for(i in 1:N_age){
    for(j in 1:N_het){
      w_HH_eq = HH_eq[i,j,]/sum(HH_eq[i,j,])
      A_par_eq_ave[i,j] = sum( A_par_eq[i,j,]*w_HH_eq )
    }
  }
  
  
  ##########################
  ## Add in maternally-acquired immunity
  
  for(j in 1:N_het){
    for(k in 1:(K_max+1)){
      # A_par_eq[,j,k] = A_par_eq[,j,k] + A_par_eq_ave[index_MI_20,j]*P_MI*exp(-age_mids/d_MI)
      A_par_mat_eq[,j,k] = A_par_eq_ave[index_MI_20,j]*P_MI*exp(-age_mids/d_MI)
    }
  }
  
  A_par_total_eq = A_par_eq + A_par_mat_eq
  
  
  ######################################################
  ## 2.3. ##  Equilibrium levels of clinical immunity ##
  ######################################################
  
  A_clin_eq  = array(NA, dim=c(N_age, N_het, K_max+1))
  A_clin_mat_eq  = array(NA, dim=c(N_age, N_het, K_max+1))
  A_clin_total_eq  = array(NA, dim=c(N_age, N_het, K_max+1))
  
  for(j in 1:N_het){
    
    #############################	
    ## Youngest age category
    G_VEC          = c(0, lam_eq[1,j]*HH_eq[1,j,1:K_max]/HH_eq[1,j,2:(K_max+1)] ) + (0:K_max)*ff
    G_VEC[K_max+1] = G_VEC[K_max+1] + lam_eq[1,j] 
    G_VEC          = G_VEC/( G_VEC*u_clin + 1 )
    
    LAM_MAT = - diag(K_max+1)
    LAM_MAT[K_max+1,K_max+1] = 0
    for(k in 1:K_max){
      LAM_MAT[k+1,k] = HH_eq[1,j,k]/HH_eq[1,j,k+1]
    }
    
    GAM_MAT = - diag(0:K_max)
    for(k in 1:K_max){
      GAM_MAT[k,k+1] = k*HH_eq[1,j,k+1]/HH_eq[1,j,k]
    }
    
    Z_MAT = - r_clin*D_MAT + lam_eq[1,j]*LAM_MAT + gamma_L*GAM_MAT - mu_H*D_MAT - r_age[1]*D_MAT
    A_clin_vec = solve( -Z_MAT )%*%( G_VEC )
    A_clin_eq[1,j,] = A_clin_vec
    
    ##########################
    ## Older age categories
    for(i in 2:N_age){
      G_VEC          = c(0, lam_eq[i,j]*HH_eq[i,j,1:K_max]/HH_eq[i,j,2:(K_max+1)] ) + (0:K_max)*ff
      G_VEC[K_max+1] = G_VEC[K_max+1] + lam_eq[i,j] 
      G_VEC          = G_VEC/( G_VEC*u_clin + 1 )
      
      r_VEC = r_age[i-1]*A_clin_eq[i-1,j,]*HH_eq[i-1,j,]/HH_eq[i,j,]
      
      LAM_MAT = - diag(K_max+1)
      LAM_MAT[K_max+1,K_max+1] = 0
      for(k in 1:K_max){
        LAM_MAT[k+1,k] = HH_eq[i,j,k]/HH_eq[i,j,k+1]
      }
      
      GAM_MAT = - diag(0:K_max)
      for(k in 1:K_max){
        GAM_MAT[k,k+1] = k*HH_eq[i,j,k+1]/HH_eq[i,j,k]
      }
      
      Z_MAT = - r_clin*D_MAT + lam_eq[i,j]*LAM_MAT + gamma_L*GAM_MAT - mu_H*D_MAT - r_age[i]*D_MAT
      A_clin_vec = solve( -Z_MAT )%*%( G_VEC + r_VEC )
      A_clin_eq[i,j,] = as.vector(A_clin_vec)
    }
  }
  
  
  ##########################
  ## Average over hypnozoites
  A_clin_eq_ave = matrix(NA, nrow=N_age, ncol=N_het)
  for(i in 1:N_age){
    for(j in 1:N_het){
      w_HH_eq = HH_eq[i,j,]/sum(HH_eq[i,j,])
      A_clin_eq_ave[i,j] = sum( A_clin_eq[i,j,]*w_HH_eq )
    }
  }
  
  ##########################
  ## Add in maternally-acquired immunity
  for(j in 1:N_het){
    for(k in 1:(K_max+1)){
      A_clin_mat_eq[,j,k] = A_clin_eq_ave[index_MI_20,j]*P_MI*exp(-age_mids/d_MI)
    }
  }
  A_clin_total_eq = A_clin_eq + A_clin_mat_eq
  
  #######################################
  ## Effects of immune functions
  r_PCR_eq  = 1/( d_PCR_min + (d_PCR_max - d_PCR_min)/( 1 + (A_par_total_eq/A_d_PCR_50pc)^K_d_PCR ) )
  phi_LM_eq = phi_LM_min + (phi_LM_max-phi_LM_min)/( 1 + (A_par_total_eq/A_LM_50pc)^K_LM ) 
  phi_D_eq  = phi_D_min  + (phi_D_max-phi_D_min)/( 1 + (A_clin_total_eq/A_D_50pc)^K_D ) 
  
  #######################################
  ## Effects of immune functions 
  ## (averaged over hypnozoites)
  phi_LM_eq_ave = matrix(NA, nrow=N_age, ncol=N_het)
  phi_D_eq_ave  = matrix(NA, nrow=N_age, ncol=N_het)
  for(i in 1:N_age){
    for(j in 1:N_het){
      w_HH_eq = HH_eq[i,j,]/sum(HH_eq[i,j,])
      phi_LM_eq_ave[i,j] = sum( phi_LM_eq[i,j,]*w_HH_eq )
      phi_D_eq_ave[i,j]  = sum( phi_D_eq[i,j,]*w_HH_eq )
    }
  }
  
  ######################################################################
  ## 2.4. ##  Function for equilibrium solution via linear algebra    ##
  ######################################################################
  
  MM_ij <- function(i, j){
    MM = matrix(0, nrow=6*(K_max+1), ncol=6*(K_max+1))
    
    MM[S_index,S_index]         = - lam_eq[i,j]*D_MAT - ff*K_MAT + gamma_L*L_MAT - mu_H*D_MAT - r_age[i]*D_MAT 
    MM[S_index,I_PCR_index]     = + t(t(D_MAT)*r_PCR_eq[i,j,])
    MM[S_index,P_index]         = + r_P*D_MAT
    
    MM[I_PCR_index,S_index]     = + lam_eq[i,j]*t(t(OD_MAT)*(1.0-phi_LM_eq[i,j,])) + ff*t(t(K_MAT)*(1.0-phi_LM_eq[i,j,]))
    MM[I_PCR_index,I_PCR_index] = - lam_eq[i,j]*D_MAT - ff*K_MAT - t(t(D_MAT)*r_PCR_eq[i,j,]) +
      lam_eq[i,j]*t(t(OD_MAT)*(1-phi_LM_eq[i,j,])) + ff*t(t(K_MAT)*(1-phi_LM_eq[i,j,])) +
      gamma_L*L_MAT - mu_H*D_MAT - r_age[i]*D_MAT
    MM[I_PCR_index,I_LM_index]  = + r_LM*D_MAT
    
    MM[I_LM_index,S_index]      = + lam_eq[i,j]*t(t(OD_MAT)*(phi_LM_eq[i,j,]*(1.0-phi_D_eq[i,j,]))) + ff*t(t(K_MAT)*(phi_LM_eq[i,j,]*(1.0-phi_D_eq[i,j,])))
    MM[I_LM_index,I_PCR_index]  = + lam_eq[i,j]*t(t(OD_MAT)*(phi_LM_eq[i,j,]*(1.0-phi_D_eq[i,j,]))) + ff*t(t(K_MAT)*(phi_LM_eq[i,j,]*(1.0-phi_D_eq[i,j,])))
    MM[I_LM_index,I_LM_index]   = - lam_eq[i,j]*D_MAT - ff*K_MAT - r_LM*D_MAT +
      lam_eq[i,j]*t(t(OD_MAT)*(1-phi_D_eq[i,j,])) + ff*t(t(K_MAT)*(1-phi_D_eq[i,j,])) +
      gamma_L*L_MAT - mu_H*D_MAT - r_age[i]*D_MAT
    MM[I_LM_index,I_D_index]    = + r_D*D_MAT
    
    MM[I_D_index,S_index]       = + lam_eq[i,j]*(1-chi_treat)*t(t(OD_MAT)*(phi_LM_eq[i,j,]*phi_D_eq[i,j,])) + ff*(1-chi_treat)*t(t(K_MAT)*(phi_LM_eq[i,j,]*phi_D_eq[i,j,]))
    MM[I_D_index,I_PCR_index]   = + lam_eq[i,j]*(1-chi_treat)*t(t(OD_MAT)*(phi_LM_eq[i,j,]*phi_D_eq[i,j,])) + ff*(1-chi_treat)*t(t(K_MAT)*(phi_LM_eq[i,j,]*phi_D_eq[i,j,]))
    MM[I_D_index,I_LM_index]    = + lam_eq[i,j]*(1-chi_treat)*t(t(OD_MAT)*phi_D_eq[i,j,]) + ff*(1-chi_treat)*t(t(K_MAT)*phi_D_eq[i,j,])
    MM[I_D_index,I_D_index]     = - lam_eq[i,j]*D_MAT - r_D*D_MAT + lam_eq[i,j]*OD_MAT +
      gamma_L*L_MAT - mu_H*D_MAT - r_age[i]*D_MAT
    
    MM[T_index,S_index]         = + lam_eq[i,j]*chi_treat*t(t(OD_MAT)*(phi_LM_eq[i,j,]*phi_D_eq[i,j,])) + ff*chi_treat*t(t(K_MAT)*(phi_LM_eq[i,j,]*phi_D_eq[i,j,]))
    MM[T_index,I_PCR_index]     = + lam_eq[i,j]*chi_treat*t(t(OD_MAT)*(phi_LM_eq[i,j,]*phi_D_eq[i,j,])) + ff*chi_treat*t(t(K_MAT)*(phi_LM_eq[i,j,]*phi_D_eq[i,j,]))
    MM[T_index,I_LM_index]      = + lam_eq[i,j]*chi_treat*t(t(OD_MAT)*phi_D_eq[i,j,]) + ff*chi_treat*t(t(K_MAT)*phi_D_eq[i,j,])
    MM[T_index,T_index]         = - lam_eq[i,j]*D_MAT - r_T*D_MAT + lam_eq[i,j]*OD_MAT + 
      gamma_L*L_MAT - mu_H*D_MAT - r_age[i]*D_MAT
    
    MM[P_index,T_index]         = + r_T*D_MAT    
    MM[P_index,P_index]         = - lam_eq[i,j]*D_MAT - r_P*D_MAT + lam_eq[i,j]*OD_MAT + 
      gamma_L*L_MAT - mu_H*D_MAT - r_age[i]*D_MAT
    
    MM
  }
  
  
  ###################################################
  ## 2.5. ##  Solve for equilibrium solution       ##
  ###################################################
  
  ###################################################
  ## Objects for storing equilibrium 
  S_eq     = array(NA, dim=c(N_age, N_het, K_max+1))
  I_PCR_eq = array(NA, dim=c(N_age, N_het, K_max+1))
  I_LM_eq  = array(NA, dim=c(N_age, N_het, K_max+1))
  I_D_eq   = array(NA, dim=c(N_age, N_het, K_max+1))
  T_eq     = array(NA, dim=c(N_age, N_het, K_max+1))
  P_eq     = array(NA, dim=c(N_age, N_het, K_max+1))
  
  
  for(j in 1:N_het){
    MM = MM_ij(1,j)
    
    BB    = rep(0, 6*(K_max+1))
    BB[1] = - w_het[j]*mu_H
    
    XX = solve(MM)%*%BB
    
    S_eq[1,j,]     = XX[S_index]
    I_PCR_eq[1,j,] = XX[I_PCR_index]
    I_LM_eq[1,j,]  = XX[I_LM_index]
    I_D_eq[1,j,]   = XX[I_D_index]
    T_eq[1,j,]     = XX[T_index]
    P_eq[1,j,]     = XX[P_index]
    
    
    for(i in 2:N_age){
      MM = MM_ij( i, j )
      BB = - r_age[i-1]*XX
      XX = solve(MM)%*%BB
      
      S_eq[i,j,]     = XX[S_index]
      I_PCR_eq[i,j,] = XX[I_PCR_index]
      I_LM_eq[i,j,]  = XX[I_LM_index]
      I_D_eq[i,j,]   = XX[I_D_index]
      T_eq[i,j,]     = XX[T_index]
      P_eq[i,j,]     = XX[P_index]
    }
  }
  
  states_3d <- list("Age" = age,
                    "HH" = HH_eq,
                    "lambda" = lam_eq,
                    "lambda_H" = lam_H_eq,
                    "S" = S_eq,
                    "U" = I_PCR_eq,
                    "A" = I_LM_eq,
                    "D" = I_D_eq,
                    "T" = T_eq,
                    "P" = P_eq,
                    "IAA" = A_par_eq,
                    "IAM" = A_par_mat_eq,
                    "ICA" = A_clin_eq,
                    "ICM" = A_clin_mat_eq,
                    "phi_clin" = phi_D_eq,
                    "phi_patent" = phi_LM_eq,
                    "r_PCR" = r_PCR_eq
  )
  
  
  #######################################################################################
  ## 3.1. ##  Vector model (copied from Nora Schmidt's deterministic equilibrium code) ##
  #######################################################################################
  
  FOIvij_eq <- array(dim=c(N_age,N_het,(K_max+1)))
  for (kk in 1:(K_max+1)){
    for (j in 1:N_het){
      for (i in 1:N_age){
        FOIvij_eq[i, j, kk] <-  aa *
          x_age_het[i,j] * (c_T * T_eq[i, j,kk] +
                              c_D *I_D_eq[i, j,kk] +
                              c_LM * I_LM_eq[i, j,kk] +
                              c_PCR * I_PCR_eq[i, j,kk])
      }
    }
  }
  
  ## When we remove the age-standardisation from the earlier EIR calculation,
  ## we need to instead account for this in the calculation of FOIM, 
  ## as found in the pf model.
  if(isTRUE(p$vivax_EIR_at_population_level)){
    FOIv_eq <- sum(FOIvij_eq)
  } else {
    FOIv_eq <- sum(FOIvij_eq) * omega_age 
  }
  
  return(list(states = states_3d, FOIM = FOIv_eq))
  
}
