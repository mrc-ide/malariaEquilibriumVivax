#------------------------------------------------
#' @title Equilibrium solution without biting heterogeneity
#'
#' @description UPDATE NOTATION: Returns the equilibrium states for the model of Griffin et al.
#'   (2014). A derivation of the equilibrium solutions can be found in Griffin
#'   (2016).
#'
#'   This function does not account for biting heterogeneity - see
#'   \code{human_equilibrium()} for function that takes this into account.
#'
#' @param EIR EIR for adults, in units of infectious bites per person per year
#' @param ft proportion of clinical cases effectively treated
#' @param p vector of model parameters
#' @param age vector of age groups, in units of years
#'
#' @references Griffin et. al. (2014). Estimates of the changing age-burden of
#'   Plasmodium falciparum malaria disease in sub-Saharan Africa.
#'   doi:10.1038/ncomms4136
#'
#'   Griffin (2016). Is a reproduction number of one a threshold for Plasmodium
#'   falciparum malaria elimination? doi:10.1186/s12936-016-1437-9 (see
#'   supplementary material)
#'
#' @export

human_equilibrium_no_het <- function(EIR, ft, p, age) {
  # check inputs
  assert_single_pos(EIR, zero_allowed = FALSE)
  assert_single_bounded(ft)
  # assert_custom_class(p, "model_params")
  assert_vector_pos(age)
  assert_noduplicates(age)
  assert_increasing(age)
  
  # get basic properties
  n_age <- length(age)
  age_days <- age*365
  EIR <- EIR/365
  p$eta <- 1/p$average_age
  p$rp <- 0.1
  
  # produce population age distribution using eta, which is defined as 1/average
  # age. The population distribution can be envisioned as an exponential
  # distribution, with mass feeding in from the left due to birth at rate eta,
  # people ageing with rates defined based on the width of the age groups, and
  # mass leaking out of all categories with death rate eta. Total birth and
  # death rates are equal, making the distribution stable.
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
      prop[i] <- p$eta/(r[i] + p$eta)
    } else {
      prop[i] <- prop[i-1]*r[i-1]/(r[i] + p$eta)
    }
  }
  
  # calculate midpoint of age range. There is no midpoint for the final age
  # group, so use beginning of range instead
  age_days_midpoint <- c((age_days[-n_age] + age_days[-1])/2, age_days[n_age])
  
  # get age category that represents a 20 year old woman
  age20 <- which.min(abs(age_days_midpoint - (20*365)))
  
  # calculate relative biting rate for each age group
  psi <- 1 - p$rho*exp(-age_days_midpoint/p$a0)
  omega <- 1 - p$rho*p$eta/(p$eta + 1/p$a0)

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

  # return matrix
  ret <- cbind(
    age = age,
    S = S,
    T = T,
    D = D,
    A = A,
    U = U,
    P = P,
    inf = inf,
    prop = prop,
    psi = psi,
    pos_LM = pos_LM,
    pos_PCR = pos_PCR,
    inc = inc,
    pat_inc = pat_inc,
    clin_inc = clin_inc,
    ICA = ICA,
    ICM = ICM,
    ID = IDA,
    IDM = IDM,
    HH = HH,
    FOI = FOI,
    FOIH = FOIH,
    phi_patent = phi_patent,
    phi_clin = phi_clin,
    r_PCR = r_PCR,
    EPS = EPS,
    r = r
  )
  return(ret)
}

#------------------------------------------------
#' @title Equilibrium solution
#'
#' @description Returns the equilibrium states for the model of Griffin et al.
#'   (2014). A derivation of the equilibrium solutions can be found in Griffin
#'   (2016). Integrates over the distribution of biting heterogeneity using
#'   Gaussian quadrature.
#'
#' @inheritParams human_equilibrium_no_het
#' @param h a list of Gauss-Hermite nodes and associated weights, used for
#'   integrating over heterogeneity in biting. See \code{?gq_normal} for an
#'   example.
#'
#' @references Griffin et. al. (2014). Estimates of the changing age-burden of
#'   Plasmodium falciparum malaria disease in sub-Saharan Africa.
#'   doi:10.1038/ncomms4136
#'
#'   Griffin (2016). Is a reproduction number of one a threshold for Plasmodium
#'   falciparum malaria elimination? doi:10.1186/s12936-016-1437-9 (see
#'   supplementary material)
#'
#' @export

human_equilibrium <- function(EIR, ft, p, age, h = gq_normal(10)) {
  
  # check inputs
  assert_single_pos(EIR, zero_allowed = FALSE)
  assert_single_bounded(ft)
  # assert_custom_class(p, "model_params")
  assert_vector_pos(age)
  assert_noduplicates(age)
  assert_increasing(age)
  assert_list(h)
  assert_in(c("nodes", "weights"), names(h))
  assert_same_length(h$nodes, h$weights)
  assert_numeric(h$nodes, h$weights)
  
  # get basic properties and initialise
  nh <- length(h$nodes)
  FOIM <- 0 		# overall force of infection on mosquitoes, weighted by onward biting rates
  
  # loop through all Gaussian quadrature nodes
  for (j in 1:nh) {
    zeta <- exp(-p$sigma_squared*0.5 + sqrt(p$sigma_squared)*h$nodes[j])
    Ej <- human_equilibrium_no_het(EIR = EIR*zeta, ft = ft, p = p, age = age)
    
    if (j == 1) {
      E <- Ej*h$weights[j]
    } else {
      E <- E + Ej*h$weights[j]
    }
    FOIM <- FOIM + sum(Ej[,"inf"]*Ej[,"psi"])*h$weights[j]*zeta
  }
  
  # calculate overall force of infection on mosquitoes
  eta <- 1/p$average_age
  omega <- 1 - p$rho*eta/(eta + 1/p$a0)
  # omega <- 1
  alpha <- p$blood_meal_rates*p$Q0
  FOIM <- FOIM*alpha/omega
  browser()
  
  # 1000*sum(E[,"inc"] * E[,"FOI"] / E[,"FOIH"])
  # 1000*sum(E[,"inc"] * (E[,"FOIH"]-E[,"FOI"]) / E[,"FOIH"])
  # 
  # 1000*sum(E[,"inc"])
  # 1000*sum(E[,"pat_inc"])
  # 1000*sum(E[,"clin_inc"])
  # 
  # sum(E[,"phi_patent"]*E[,"prop"])
  # sum(E[,"phi_clin"]*E[,"prop"])
  # 
  # 1000*sum(E[,"HH"]*E[,"prop"])*p$f
  
  round(colSums(E[,c("S","D","A","U","T","P")]),3)
  
  round(sum(E[,"HH"]*E[,"prop"]),3)
  round(sum(E[,"ICA"]*E[,"prop"]),1)
  round(sum(E[,"ID"]*E[,"prop"]),1)
  
  
  # return as list
  return(list(
    states = E,
    FOIM = FOIM
  ))
}
