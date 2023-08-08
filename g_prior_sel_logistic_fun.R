# G Prior function --------------------------------
g.prior.sel.logistic = function(X=NULL, Y=NULL, U=NULL, SETNUM = NULL){
  # G Prior Model ----------------------------------------------
  g_prior_sel.model <- 
    "model {

  for(i in 1:N) {
    Y[i] ~ dbern(mu[i])
    logit(mu[i]) <- alpha + inprod(b[1:P], X[i,1:P]) + inprod(delta[1:Q], U[i,1:Q])
  }

  # prior on covariate effects
  for(q in 1:Q) { delta[q] ~ dnorm(0, 1.0E-06) }
  
  # prior on intercept
  alpha ~ dnorm(0, 1.0E-06)
  
  # prior on exposure effects
  beta[1:P] ~ dmnorm(mu.beta[1:P], T[1:P, 1:P])
  for(j in 1:P) {
    mu.beta[j] <- (1-gamma[j])*prop.mu.beta[j]
    b[j] <- beta[j]*gamma[j]
    gamma[j] ~ dbern(pi)
    for(k in 1:P) {
      T[j,k] <- gamma[j]*gamma[k]*XtX[j,k]/(G) + (1-gamma[j]*gamma[k])*equals(j,k)*pow(prop.sd.beta[j],-2)
    }
  }
  pi ~ dbeta(1,P)
  #pi ~ dbeta(P, 1)

  # semi-Bayes
  # G <- w/(1-w)
  # w <- .99   # w -> 0 shrink to common mean; as w -> inf toward the maximum likelihood estimate

  # Zellner and Siow prior on G
  #b0 <- 0.5*N
  #inv.G ~ dgamma(0.5, b0)
  #G <- 1/inv.G
  #w <- G/(G+1)

  # Hyper-g prior (following Perrakis 2018, note that this is on the G^-1 so the Beta distribution is switchd in terms of a and b from Li and Clyde 2019 equation 34)
  a <- 3
  bw <- a/2 - 1
  w~dbeta(1,bw)
  G <- w/(1-w)

  # Hyper-g/n prior (following Perrakis 2018, note that this is on the G^-1 so the Beta distribution is switchd in terms of a and b from Li and Clyde 2019 equation 34)
  #a <- 3
  #bw <- a/2 - 1
  #w~dbeta(1,bw)
  #G <- N*w/(1-w)

  #beta-prime 
  #G <- w/(1-w)
  #w ~ dbeta(bw, .25)
  #bw <- (N-P_m-1.5)/2
  #P_m <- sum(gamma[1:P])

  # g-estimation
  eta.low <- inprod(b[1:P], profiles[1,1:P])
  eta.high <- inprod(b[1:P], profiles[2,1:P])
  psi <-eta.high-eta.low
  
  }"
  
  # Other Stuff -------------------------------------------------------------
  N <- length(Y)
  P <- ncol(X)
  Q <- ncol(U)
  exposure.Names = colnames(X)
  # Set profiles
  profiles = c(-1,1)*matrix(.5, nrow=2, ncol=P)
  exposure.Names = colnames(X)
  
  
  ### get the univariate result
  univariate.results <- t(
    sapply(1:P, 
           FUN=function(p) {  # using index p facilitate write
             x <- as.matrix(X[,p])
             reg <- clogit(Y~x + strata(SETNUM))      # perform logistic regression
             s.reg <- summary(reg)                 # get the summary for the regression
             c.reg <- s.reg$coef["x",]             # select the coefficients for the exposure
             return(c.reg)                         
           }, 
           simplify=T))
  
  univariate.results <- data.frame(exposure.Names,
                                   univariate.results)
  
  ### g prior model result
  prop.mu.beta <- univariate.results$coef
  prop.sd.beta <- univariate.results$se.coef.
  XtX <- t(as.matrix(X))%*%as.matrix(X) 
  # Data for Jags Model
  jdata <- list(N=N, Y=Y, X=X, U=U, P=P, Q=Q, XtX=XtX, 
                profiles=profiles,
                prop.mu.beta=prop.mu.beta, 
                prop.sd.beta=prop.sd.beta)
  # Variables to sample
  var.s <- c("b", "beta","gamma", "psi", "pi", "G", "w")
  # Fit Model
  model.fit <- jags.model(file=textConnection(g_prior_sel.model), 
                          data=jdata, n.chains=1, n.adapt=4000, quiet=T)
  # Update Model
  update(model.fit, n.iter=1000, progress.bar="none")
  model.fit <- coda.samples(model=model.fit,
                            variable.names=var.s, 
                            n.iter=50000, 
                            thin=1, 
                            progress.bar="none")
  
  r <- summary(model.fit)
  var.names <- c("G", 
                 paste(exposure.Names, "b", sep="."),
                 paste(exposure.Names, "beta", sep="."), 
                 paste(exposure.Names, "gamma", sep="."),
                 "pi",
                 "psi",
                 "w")
  g_prior.results <- data.frame(var.names, 
                                r$statistics[,1:2], 
                                r$quantiles[,c(1,5)])
  
  # Calculate Wald statistic
  g_prior.results$wald = ifelse(str_detect(g_prior.results$var.names, 
                                           'beta') | 
                                  str_detect(g_prior.results$var.names, 
                                             'psi'), 
                                abs(g_prior.results[,"Mean"]/
                                      g_prior.results[,"SD"]), 
                                NA)
  # Calculate p value
  g_prior.results$p.val = (2*(1-pnorm(g_prior.results$wald, 0, 1)))
  
  return(g_prior.results)
}