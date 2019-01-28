
# Fit BMA coefficient for a beta distribution with discrete component at 1
#' @param train.tel Vector of training telemetry data on [0,1]
#' @param train.ens Matrix of training ensemble member data [time x member] on [0,1]
#' @return A formula for a discrete-continuous mixture model with beta distribution
beta1_model <- function(tel, ens) {

  # 1. Logistic regression for a's
  # Define what the the observations should be here
  discrete_model <- summary(glm(obs ~ fc,family=binomial(link='logit'),data=data.frame(fc= , obs = )))

  # 2. linear regression for b's
  # Define what the the observations should be here
  mean_model <- summary(lm(obs ~ fc, data=data.frame(fc= , obs= )))
  # 3. ME algorithm for w's and c's

  # 4. back out p, y, and define final discrete-continuous model
 }


# Expectation-maximization function, modified from code courtesy of Will Kleiber
em <- function(FCST, OBS, A0, A1, A2, B0, B1, C0=0.5, C1=0.5, eps=1e-005, maxiter=1000, start.w=NULL)

    # MODEL for one forecast : y_s is solar power, logit P(y_s = 1 | f) = a0 + a1 f
    #                          solar power level, conditional on it being less than rated power (i.e., 1) is beta distributed
    #                          with mean = b0 + b1 f, standard deviation = c0 + c1 f

    # Inputs:
    #  FCST        array of dimension (length of training period)x(number of stations)x(number of ensemble members)
    #              where each matrix (FCST[,,k]) has a ensemble member's forecast for these sites and times
#  OBS         matrix of dimension (length of training period)x(number of stations)
#              of observations at these sites and times
#  XN          for X=A,B and N=0,1,2, array of dimension (length of training period)x(number of stations)x
#              (number of ensemble members) of estimates of XN (usually posterior mean for A0,1,2 and
#              least squares from regression for B0,1)
#               Can set A2's to 0's to ignore the indicator function aspect
#  CN          for N=0,1 starting estimates of CN, which are assumed constant across sites and members (equal variances among ensemble member)

#  eps     stopping criterion
#  maxiter maximum number of EM iterations allowed
#  start.w initial values for the weights (optional)

{
  # number of ensemble members
  K <- dim(FCST)[3]
  ntime <- dim(FCST)[1]
  nsite <- dim(FCST)[2]

  # set intial weights, either as equal weights or given starting values, ignoring forecast members which are unavailable
  avail <- apply(FCST, MARGIN = 3, FUN = function(x) {!all(is.na(x))}) # Which members are available
  if(is.null(start.w)){
    w <- sapply(avail, function(x) ifelse(x, 1/sum(avail), 0))
  }
  else{
    w <- sapply(seq_along(start.w), function(k) ifelse(avail[k], start.w[k]/sum(start.w[avail]), 0))
  }

  error <- rep(1,times=(K))
  count <- 0

  # Precalculate probability of clipping and beta density
  PoC <- sapply(1:ntime, function(t) {sapply(1:nsite, function(s) {sapply(1:K, get_poc(FCST[t,s,k], A0[t,s,k], A1[t,s,k], A2[t,s,k]))}, simplify="matrix")}, simplify="array")


  # main EM algorithm
  while((max(abs(error)) > eps) && (count < maxiter))
  {
    # Re-calculate beta density estimates based on current estimate for C0 and C1
    db <- sapply(1:ntime, function(t) {sapply(1:nsite, function(s) {sapply(1:K, get_db(OBS[t,s,k], FCST[t,s,k], B0[t,s,k], B1[t,s,k], C0, C1))}, simplify="matrix")}, simplify="array")

    ## E step
    # sumz is weighted sum of density functions for each member
    # PoC is a vector of probabilities of clipping, length of the training data

    # z is the E step of algorithm, estimated for each site and time and forecast
    # z is an array with the first entry being day, second entry site, third entry forecast
    z_num <- sapply(1:nsite, function(s) {
      sapply(1:K, function(k) {get_z_num_by_member(OBS[,s], FCST[,s,k], PoC[,s,k], db[,s,k], w[k])},
             simplify="matrix")},simplify='array')

    # TODO Check dimensions -- aperm??
    # TODO TO DECIDE: SHOULD sites BE A SEPARATE MARGIN (allowing local A's and B's and global C's) OR NOT? Local A'S AND B'S CAN ALSO BE REPLICATED DOWN THE TRAINING DATA VECTOR TO ACCOMPLISH THE SAME THING
    sumz <- apply(z_num, MARGIN=2, FUN=sum, na.rm=T)
    z <- sapply(1:nsite, function(s) {sapply(1:K, function(k) {z_num[s,k,]/sumz[k]}, simplify="matrix")}, simplify="array" )

    ## M step
    # new weights and variance deflation
    w.new <- sapply(1:K, function(k) {ifelse(any(!is.na(z[,,k])), mean(z[,,k], na.rm=TRUE), 0)})

    # change from previous iteration to this one
    error <- sapply(1:K, function(k) {return(w.new[k] - w[k])})

    w <- w.new

    # Using the new w estimate, re-optimize C0 and C1
    get_log_lik <- function(params) {
      # THESE ARE REDUNDANT WITH ABOVE
      db <- sapply(1:ntime, function(t) {sapply(1:nsite, function(s) {sapply(1:K, get_db(OBS[t,s,k], FCST[t,s,k], B0[t,s,k], B1[t,s,k], params[1], params[2]))}, simplify="matrix")}, simplify="array")
      z_num <- sapply(1:nsite, function(s) {
        sapply(1:K, function(k) {get_z_num_by_member(OBS[,s], FCST[,s,k], PoC[,s,k], db[,s,k], w[k])},
               simplify="matrix")},simplify='array')
      sumz <- apply(z_num, MARGIN=2, FUN=sum, na.rm=T)
      return(sum(log(sumz), na.rm=T))
    }
    optim_list <- optim(c(C0, C1), get_log_lik, control = list(fnscale=-1))
    if (optim_list$convergence != 0) stop(optim_list$message)

    error <- append(error, optim_list$par - c(C0, C1)) # Complete list of changes to all w and c parameters
    C0 <- optim_list$par[1]
    C1 <- optim_list$par[2]

    count <- count + 1
  }

  # Compute log-likelihood (corrected by Adrian on 10/21/03)
  lik=sum(log(sumz), na.rm=T)
  return(list(loglik=lik, w=w, count=count))
}



# Get time-series of z numerator for each member
# Returns a time-series vector of z for a given member at a given site
get_z_num_by_member <- function(OBS, FCST, PoC, db, w){
  z <- rep(NA, length(FCST))
  notna <- !is.na(OBS) & !is.na(FCST)
  # Check whether this ensemble member is available
  if (any(notna)){
    z[notna] <- mapply(FUN=get_z, OBS=OBS[notna], PoC=PoC[notna], db=db[notna], w=w)
  }
  return(z)
}

# Get z for single instance
# density is (PoC)*1[obs==1] + (1-PoC)*Beta(obs,a,b)*1[obs < 1]
get_z <- function(OBS, PoC, db, w) {
  ifelse(OBS==1, w*PoC, w*(1-PoC)*db)
}

# Get PoC for single instance
# PoC is estimated from logit(PoC)=a0 + a1*f + a2*1[obs==1]
# This seems inverted from the PoP defintion from Will Kleiber's code, which reads to me like it should be 1-PoP
get_poc <- function(FCST, A0, A1, A2) {
  if (is.na(FCST)) {return(NA)}
  else return(exp(A0+A1*FCST+A2*(FCST==1))/(1+exp(A0+A1*FCST+A2*(FCST==1))))
}

# Get beta density for single instance
get_beta_density <- function(OBS, FCST, B0, B1, C0, C1) {
  if (is.na(FCST) | is.na(OBS)) {return(NA)}
  else if (OBS==1) {return(NA)}
  else {
    # Estimate mu (mean)
    mu <- B0 + B1*FCST
    # Estimate sigma (std deviation)
    sigma <- C0 + C1*FCST
    # Back-calculate gamma (I'm not sure this works out to be positive, given my current algebra) *****
    gamma <- mu*(1-mu)/sigma^2 - 1
    if (gamma <=0) stop(paste("Non-positive gamma in beta density estimation with mu", mu, "and sigma", sigma, sep=" "))

    # Back-calculate alpha and beta from gamma and rho=mu
    alpha <- mu*gamma
    beta <- gamma *(1-mu)

    return(dbeta(OBS,shape1=alpha,shape2=beta))
  }
}
