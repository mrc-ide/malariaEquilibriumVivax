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

vivax_equilibrium <- function(age, ft, EIR, p, v_eq = "full"){
  
  ## Check Parameters
  if(!is.numeric(age)) stop("age provided is not numeric")
  if(!is.numeric(p$n_heterogeneity_groups)) stop("het_brackets provided is not numeric")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  if(!is.numeric(EIR)) stop("EIR provided is not numeric")
  if(!is.numeric(p$kmax)) stop("K_max is not provided or not numeric")
  if(!(p$kmax %in% c(1:30))) stop("K_max is not an integer between 1 and 30")
  if(!v_eq %in% c("full","simplified")){stop("vivax equilibrium must be 'full' or 'simplified'")}
  
  assert_single_pos(EIR, zero_allowed = FALSE)
  assert_single_bounded(ft)
  assert_custom_class(p, "list")
  assert_vector_pos(age)
  assert_noduplicates(age)
  assert_increasing(age)
  
  if(v_eq == "simplified"){
    return(vivax_equilibrium_simplified(age, ft, EIR, p))
  } else {
    
    
    library(survival)
    library(statmod)
    library(binom)
    
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
    
    age_bite = 1 - rho_age*exp( -age_mids/age_0 )
    omega_age = 1/sum( age_demog*age_bite )
    age_bite = omega_age*age_bite
    
    ###########################################################################
    ## Heterogeneity in mosquito bites
    
    sig_het = sqrt(sigma_squared)
    x_het <- exp(gauss.quad.prob(N_het, dist="normal", mu=-0.5*sig_het^2, sigma=sig_het)$nodes)
    w_het <-     gauss.quad.prob(N_het, dist="normal", mu=-0.5*sig_het^2, sigma=sig_het)$weights
    
    ##########################################################
    ## Not sure if GQ works perfectly with log-Normal
    ## for small N_het - might need to switch to Gamma
    ## x_het <- gauss.quad.prob(N_het, dist="gamma", alpha=1, beta=1)$nodes
    ## w_het <- gauss.quad.prob(N_het, dist="gamma", alpha=1, beta=1)$weights
    
    x_age_het <- age_bite%o%x_het
    w_age_het <- age_demog%o%w_het
    
    ################################################################
    ## 1.2. ##  Parameter definitions and additional assignments  ##
    ################################################################
    
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
    
    #########################################################
    ## 2.6. ##  Age stratify                               ##
    #########################################################
    
    
    lam_mosq <- c_PCR*sum(I_PCR_eq) + 
      c_LM*sum(I_LM_eq) + 
      c_D*sum(I_D_eq) + 
      c_T*sum(T_eq) 
    
    E_M = (lam_mosq/(mu_M + lam_mosq))*(1 - exp(-mu_M*tau_M))
    
    I_M = (lam_mosq/(mu_M + lam_mosq))*exp(-mu_M*tau_M)
    
    inf_rvoir = rep(NA, N_age)
    
    for(i in 1:N_age)
    {
      inf_rvoir[i] = c_PCR*sum(I_PCR_eq[i,,]*HH_eq[i,,]) + 
        c_LM*sum(I_LM_eq[i,,]*HH_eq[i,,]) + 
        c_D*sum(I_D_eq[i,,]*HH_eq[i,,]) + 
        c_T*sum(T_eq[i,,]*HH_eq[i,,])
    }
    
    inf_rvoir <- inf_rvoir*(age_bite/age_demog)/sum(inf_rvoir*(age_bite/age_demog))
    
    av0 <- Q0 * 1 / (foraging_time + 1/blood_meal_rates)
    mv0 <- sum(omega_age * age_demog) * EIR_site/(I_M * av0)
    
    states <- lapply(1:N_het, function(het){
      data.frame("age" = age_bounds[-N_age-1]/365,
                 "prop" = age_demog,
                 "S" = rowSums(S_eq[,het,]),
                 "U" = rowSums(I_PCR_eq[,het,]),
                 "A" = rowSums(I_LM_eq[,het,]),
                 "D" = rowSums(I_D_eq[,het,]),
                 "T" = rowSums(T_eq[,het,]),
                 "P" = rowSums(P_eq[,het,]),
                 "ID" = rowSums(A_par_eq[,het,]*HH_eq[,het,])/rowSums(HH_eq[,het,]),
                 "IDM" = rowSums(A_par_mat_eq[,het,]*HH_eq[,het,])/rowSums(HH_eq[,het,]),
                 "ICA" = rowSums(A_clin_eq[,het,]*HH_eq[,het,])/rowSums(HH_eq[,het,]),
                 "ICM" = rowSums(A_clin_mat_eq[,het,]*HH_eq[,het,])/rowSums(HH_eq[,het,]),
                 "HH" = rowSums(t(t(HH_eq[,het,]/rowSums(HH_eq[,het,]))*0:10)),
                 "EIR" = EIR_site,
                 "inf" = lam_mosq,
                 "n_hypnozoites" = 1 - HH_eq[,het,1]/rowSums(HH_eq[,het,]),
                 "E_M" = E_M,
                 "I_M" = I_M,
                 "inf_rvoir" = inf_rvoir,
                 "mv0" = mv0,
                 "psi" = age_bite,
                 "phi_clin" = rowSums(phi_D_eq[,het,]*HH_eq[,het,])/rowSums(HH_eq[,het,]),
                 "phi_patent" = rowSums(phi_LM_eq[,het,]*HH_eq[,het,])/rowSums(HH_eq[,het,]),
                 "rPCR" = rowSums(r_PCR_eq[,het,]*HH_eq[,het,])/rowSums(HH_eq[,het,]))
    })
    
    #######################################################################################
    ## 3.1. ##  Vector model (copied from Nora Schmidt's deterministic equilibrium code) ##
    #######################################################################################
    
    FOIvij_eq <- array(dim=c(N_age,N_het,(K_max+1)))
    for (kk in 1:(K_max+1)){
      for (j in 1:N_het){
        for (i in 1:N_age){
          FOIvij_eq[i, j, kk] <-  blood_meal_rates * Q0 * 
            x_age_het[i,j] * (c_T * T_eq[i, j,kk] +
                                c_D *I_D_eq[i, j,kk] +
                                c_LM * I_LM_eq[i, j,kk] +
                                c_PCR * I_PCR_eq[i, j,kk])
        }
      }
    }
    
    # Mosquito states
    FOIv_eq <- sum(FOIvij_eq)
    
    detach(p)
    return(list(states = states, FOIM = FOIv_eq))
    
  }
}

#------------------------------------------------
#' @title Vivax approximate equilibrium solution with biting heterogeneity
#'
#' @description Returns the simplified vivax equilibrium states for the model of White et al.,
#'   (2018). This solution does not loop over hypnozoite batches, instead calculating an average
#'   for each age and biting group.
#'
#' @param EIR EIR for adults, in units of infectious bites per person per year
#' @param ft proportion of clinical cases effectively treated
#' @param p vector of model parameters
#' @param age vector of age groups, in units of years
#'
#' @export

vivax_equilibrium_simplified <- function(age, ft, EIR, p, v_eq = "full"){
  
  library(survival)
  library(statmod)
  library(binom)
  
  suppressMessages(attach(p)) 
  
  #####################################
  ## 2.1. ##  Age and heterogeneity  ##
  #####################################  
  
  ##################################### 
  ## Age demography
  
  max_age  = 80*365
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
  
  for(i in 1:(N_age-1))
  {
    r_age[i] <- ( 1 - sum(age_demog[1:i]) )/( age_demog[i]*mean_age )
  }	
  
  r_age[N_age] = 0
  
  
  ###################################
  ##
  ## Modification of aging rates for immune functions.
  ##
  ## Proportions infected must always accord to the demographic
  ## constraints when stratified. Immunity is a scalar and can 
  ##(in theory) increase without bound. We must make an adjustment
  ## due to demography for aging of immunity.
  ## Provide better explanation.
  
  r_imm_age <- rep(NA, N_age)
  
  r_imm_age[1:(N_age-1)] <- r_age[1:(N_age-1)]*(age_demog[1:(N_age-1)]/age_demog[2:N_age])
  r_imm_age[N_age] <- 0
  
  
  
  ###########################################################################
  ## Age-dependent exposure to mosquito bites
  
  age_bite = 1 - rho_age*exp( -age_mids/age_0 )
  
  omega_age = 1/sum( age_demog*age_bite )
  
  age_bite = omega_age*age_bite
  
  
  ###########################################################################
  ## Heterogeneity in mosquito bites
  
  sig_het = sqrt(sigma_squared)
  
  x_het <- exp(gauss.quad.prob(N_het, dist="normal", mu=-0.5*sig_het^2, sigma=sig_het)$nodes)
  w_het <-     gauss.quad.prob(N_het, dist="normal", mu=-0.5*sig_het^2, sigma=sig_het)$weights
  
  
  ##########################################################
  ## Not sure if GQ works perfectly with log-Normal
  ## for small N_het - might need to switch to Gamma
  ## x_het <- gauss.quad.prob(N_het, dist="gamma", alpha=1, beta=1)$nodes
  ## w_het <- gauss.quad.prob(N_het, dist="gamma", alpha=1, beta=1)$weights
  
  
  x_age_het <- age_bite%o%x_het
  
  w_age_het <- age_demog%o%w_het
  
  
  
  #######################################
  ## 2.2. ##  Define model parameters  ##
  #######################################  
  
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
    
    
    ###################################################
    ## 3.2. ##  FoI and proportion with hypnozoites  ##
    ###################################################
    
    lam_eq = EIR_site*bb*x_age_het
    
    HH_eq <- matrix(NA, nrow=N_age, ncol=N_het)
    
    for(j in 1:N_het)
    {
      HH_eq[1,j] = ( 1/(gamma_L + r_imm_age[1]) )*( lam_eq[1,j] )
      
      for(i in 2:N_age)
      {
        HH_eq[i,j] = ( 1/(gamma_L + r_imm_age[i]) )*( lam_eq[i,j] + r_imm_age[i]*HH_eq[i-1,j] )
      }
    }
    
    lam_H_eq = lam_eq + ff*HH_eq
    
    
    ###################################################
    ## 3.3. ##  Equilibrium levels of immunity       ##
    ###################################################
    
    A_par_eq <- matrix(NA, nrow=N_age, ncol=N_het)
    A_clin_eq <- matrix(NA, nrow=N_age, ncol=N_het)
    
    for(j in 1:N_het)
    {
      A_par_eq[1,j] = ( 1/(r_par + r_imm_age[1]) )*( lam_H_eq[1,j]/(lam_H_eq[1,j]*u_par+1) )
      
      for(i in 2:N_age)
      {
        A_par_eq[i,j] = ( 1/(r_par + r_imm_age[i]) )*( lam_H_eq[i,j]/(lam_H_eq[i,j]*u_par+1) + r_imm_age[i]*A_par_eq[i-1,j] )
      }
    }
    
    
    for(j in 1:N_het)
    {
      A_clin_eq[1,j] = ( 1/(r_clin + r_imm_age[1]) )*( lam_H_eq[1,j]/(lam_H_eq[1,j]*u_clin+1) )
      
      for(i in 2:N_age)
      {
        A_clin_eq[i,j] = ( 1/(r_clin + r_imm_age[i]) )*( lam_H_eq[i,j]/(lam_H_eq[i,j]*u_clin+1) + r_imm_age[i]*A_clin_eq[i-1,j] )
      }
    }
    
    A_par_mat_eq <- ( P_MI*exp(-age_mids/d_MI) )%o%A_par_eq[index_MI_20,]
    A_clin_mat_eq <- ( P_MI*exp(-age_mids/d_MI) )%o%A_clin_eq[index_MI_20,]
    
    A_par_total_eq  = A_par_eq + A_par_mat_eq
    A_clin_total_eq = A_clin_eq + A_clin_mat_eq
    
    
    #######################################
    ## Effects of immune functions
    
    r_PCR_eq  = 1/( d_PCR_min + (d_PCR_max - d_PCR_min)/( 1 + (A_par_total_eq/A_d_PCR_50pc)^K_d_PCR ) )
    phi_LM_eq = phi_LM_min + (phi_LM_max-phi_LM_min)/( 1 + (A_par_total_eq/A_LM_50pc)^K_LM ) 
    phi_D_eq  = phi_D_min  + (phi_D_max-phi_D_min)/( 1 + (A_clin_total_eq/A_D_50pc)^K_D ) 
    
    
    
    ###################################################
    ## 3.4. ##  Function for equilibrium solution    ##
    ##      ##  via linear algebra                   ##
    ###################################################
    
    MM_ij <- function(i, j)
    {
      MM <- matrix(0, nrow=6, ncol=6)
      
      MM[1,1] = - lam_H_eq[i,j] - mu_H - r_age[i]
      MM[1,2] = + r_PCR_eq[i,j]
      MM[1,6] = + r_P
      
      MM[2,1] = + lam_H_eq[i,j]*(1.0-phi_LM_eq[i,j])
      MM[2,2] = + lam_H_eq[i,j]*(1.0-phi_LM_eq[i,j]) - lam_H_eq[i,j] - r_PCR_eq[i,j] - mu_H - r_age[i]	
      MM[2,3] = + r_LM
      
      MM[3,1] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*(1.0-phi_D_eq[i,j])
      MM[3,2] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*(1.0-phi_D_eq[i,j])
      MM[3,3] = + lam_H_eq[i,j]*(1.0-phi_D_eq[i,j]) - lam_H_eq[i,j] - r_LM - mu_H - r_age[i]
      MM[3,4] = + r_D
      
      MM[4,1] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*phi_D_eq[i,j]*(1-chi_treat)
      MM[4,2] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*phi_D_eq[i,j]*(1-chi_treat)
      MM[4,3] = + lam_H_eq[i,j]*phi_D_eq[i,j]*(1-chi_treat)
      MM[4,4] = - r_D - mu_H - r_age[i]
      
      MM[5,1] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*phi_D_eq[i,j]*chi_treat
      MM[5,2] = + lam_H_eq[i,j]*phi_LM_eq[i,j]*phi_D_eq[i,j]*chi_treat
      MM[5,3] = + lam_H_eq[i,j]*phi_D_eq[i,j]*chi_treat
      MM[5,5] = - r_T - mu_H - r_age[i]
      
      MM[6,5] = + r_T
      MM[6,6] = - r_P - mu_H - r_age[i]
      
      MM
    }
    
    
    ###################################################
    ## 3.5. ##  Solve for equilibrium solution       ##
    ###################################################
    
    ###################################################
    ## Objects for storing equilibrium 
    
    S_eq     <- matrix(NA, nrow=N_age, ncol=N_het)
    I_PCR_eq <- matrix(NA, nrow=N_age, ncol=N_het)
    I_LM_eq  <- matrix(NA, nrow=N_age, ncol=N_het)
    D_eq     <- matrix(NA, nrow=N_age, ncol=N_het)
    T_eq     <- matrix(NA, nrow=N_age, ncol=N_het)
    P_eq     <- matrix(NA, nrow=N_age, ncol=N_het)
    
    
    for(j in 1:N_het)
    {
      MM <- MM_ij(1,j)
      
      BB    = rep(0, 6)
      BB[1] = - w_het[j]*mu_H
      
      XX = solve(MM)%*%BB
      
      S_eq[1,j]     = XX[1]
      I_PCR_eq[1,j] = XX[2]
      I_LM_eq[1,j]  = XX[3]
      D_eq[1,j]     = XX[4]
      T_eq[1,j]     = XX[5]
      P_eq[1,j]     = XX[6]
      
      
      for(i in 2:N_age)
      {
        MM <- MM_ij( i, j )
        
        BB = - r_age[i-1]*XX
        
        XX = solve(MM)%*%BB
        
        S_eq[i,j]     = XX[1]
        I_PCR_eq[i,j] = XX[2]
        I_LM_eq[i,j]  = XX[3]
        D_eq[i,j]     = XX[4]
        T_eq[i,j]     = XX[5]
        P_eq[i,j]     = XX[6]
      }
    }
    
    

    #########################################################
    ## 3.7. ##  Mosquitoes and infectious reservoirs       ##
    #########################################################
    
    lam_mosq <- c_PCR*sum(I_PCR_eq) + 
      c_LM*sum(I_LM_eq) + 
      c_D*sum(D_eq) + 
      c_T*sum(T_eq) 
    
    E_M = (lam_mosq/(mu_M + lam_mosq))*(1 - exp(-mu_M*tau_M))
    
    I_M = (lam_mosq/(mu_M + lam_mosq))*exp(-mu_M*tau_M)
    
    
    inf_rvoir <- c_PCR*rowSums(t(t(I_PCR_eq)*x_het)) + 
      c_LM*rowSums(t(t(I_LM_eq)*x_het)) + 
      c_D*rowSums(t(t(D_eq)*x_het)) + 
      c_T*rowSums(t(t(T_eq)*x_het)) 
    
    inf_rvoir <- inf_rvoir*(age_bite/age_demog)/sum(inf_rvoir*(age_bite/age_demog))
    
    ##################################
    ##### Add FOIM calculations ######
    ##################################
    
    #########################################################
    ##  3.8. ##  Summary output                            ##
    #########################################################
    
    #########################################################
    ##  Return output
    
    inf_rvoir <- inf_rvoir*(age_bite/age_demog)/sum(inf_rvoir*(age_bite/age_demog))
    
    av0 <- Q0 * 1 / (foraging_time + 1/blood_meal_rates)
    mv0 <- sum(omega_age * age_demog) * EIR_site/(I_M * av0)
    
    states <- lapply(1:N_het, function(het){
      data.frame("age" = age_bounds[-N_age-1]/365,
                 "prop" = age_demog,
                 "S" = S_eq[,het],
                 "U" = I_PCR_eq[,het],
                 "A" = I_LM_eq[,het],
                 "D" = D_eq[,het],
                 "T" = T_eq[,het],
                 "P" = P_eq[,het],
                 "ID" = A_par_eq[,het],
                 "IDM" = A_par_mat_eq[,het],
                 "ICA" = A_clin_eq[,het],
                 "ICM" = A_clin_mat_eq[,het],
                 "HH" = HH_eq[,het],
                 "EIR" = EIR_site,
                 "inf" = lam_mosq,
                 # "n_hypnozoites" = 1 - HH_eq[HH_eq<=1] <- 1,
                 "E_M" = E_M,
                 "I_M" = I_M,
                 "inf_rvoir" = inf_rvoir,
                 "mv0" = mv0,
                 "psi" = age_bite,
                 "phi_clin" = phi_D_eq[,het],
                 "phi_patent" = phi_LM_eq[,het],
                 "dPCR" = r_PCR_eq[,het])
    })
    
    #######################################################################################
    ## 3.1. ##  Vector model (adapted from Nora Schmidt's deterministic equilibrium code) ##
    #######################################################################################
    
    FOIvij_eq <- array(dim=c(N_age,N_het))
      for (j in 1:N_het){
        for (i in 1:N_age){
          FOIvij_eq[i, j] <-  blood_meal_rates * Q0 * 
            x_age_het[i,j] * (c_T * T_eq[i, j] +
                                c_D * D_eq[i, j] +
                                c_LM * I_LM_eq[i, j] +
                                c_PCR * I_PCR_eq[i, j])
        }
      }
    
    # Mosquito states
    FOIv_eq <- sum(FOIvij_eq)

    detach(p)
    return(list(states = states, FOIM = FOIv_eq))
    
}
