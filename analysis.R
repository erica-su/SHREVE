# load packages
library(parallel)
library(tidyverse)
library(nimble)
library(rstudioapi)

# SHREVE model for simulated AGPS data

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

#---- Read Data ----
SP_data = readRDS("data.rds") # simulated AGPS data
init_data = readRDS("init_data.rds") # initial values

# set up parallel cluster (9 cores)
this_cluster <- makeCluster(9, setup_strategy = "sequential")

# ---- Nimble Function ----
run_MCMC_allcode <- function(seed, SP_data, init_data) {
  # load packages
  library(nimble)
  library(tidyverse)
  
  # variables for model
  n = dim(SP_data)[1] # number of observations
  n_id = n_distinct(SP_data$id_num) # number of subjects
  id = SP_data$id_num # ID numbers
  fuyrs = SP_data$fuyrs # follow-up time
  y = SP_data$GCC # GCC thickness
  sp = SP_data$sp_num # superpixel number
  n_sp = n_distinct(sp) # number of superpixels
  dist_data = SP_data %>% distinct(sp_num, x, y) # x-y coordinates for distances
  dists = as.matrix(dist(dist_data[, c("x", "y")], diag = T, upper = T))
  sp2 = sp + n_sp # superpixel indexing for slopes
  sp3 = sp2 + n_sp # superpixel indexing for log residual SDs
  
  n_i = SP_data %>% distinct(id_num, fuyrs) %>%
    group_by(id_num) %>%
    mutate(visit = row_number()) %>% ungroup() %>%
    mutate(visit_num = row_number())
  n_visits = dim(n_i)[1] # number of visits
  visit_num = (SP_data %>% 
                 left_join(n_i, by = c("id_num", "fuyrs")))$visit_num
  
  # exponential covariance kernel function for GP
  expcov <- nimbleFunction(
    run = function(dists = double(2), rho = double(0), sigma = double(0)) {
      returnType(double(2))
      n <- dim(dists)[1]
      result <- matrix(nrow = n, ncol = n, init = FALSE)
      sigma2 <- sigma * sigma
      for(i in 1:n)
        for(j in 1:n)
          result[i, j] <- sigma2 * exp(-dists[i, j] / rho)
      return(result)
    })
  cExpcov <- compileNimble(expcov)
  
  # exponential covariance kernel function for MGP
  covCombine <- nimbleFunction(
    run = function(dists = double(2), rho1 = double(0), rho2 = double(0),
                   rho3 = double(0), rho12 = double(0), rho13 = double(0),
                   rho23 = double(0), sigma1 = double(0), sigma2 = double(0), 
                   sigma3 = double(0), R_v12 = double(0), R_v13 = double(0),
                   R_v23 = double(0)) {
      returnType(double(2))
      n <- dim(dists)[1]
      a <- n + 1
      b <- 2 * n
      c <- b + 1
      d <- 3 * n
      sigma1_sq <- sigma1 * sigma1
      sigma2_sq <- sigma2 * sigma2
      sigma3_sq <- sigma3 * sigma3
      rss12 <- R_v12 * sigma1 * sigma2 * rho12 / sqrt(rho1 * rho2)
      rss13 <- R_v13 * sigma1 * sigma3 * rho13 / sqrt(rho1 * rho3)
      rss23 <- R_v23 * sigma2 * sigma3 * rho23 / sqrt(rho2 * rho3)
      result <- matrix(nrow = d, ncol = d, init = FALSE)
      for(i in 1:n)
        for(j in 1:n)
          result[i, j] <- sigma1_sq * exp(-dists[i, j] / rho1)
      for(i in 1:n)
        for(j in a:b) {
          result[i, j] <- rss12 * exp(-dists[i, j - n] / rho12)
          result[j, i] <- result[i, j]
        }
      for(i in 1:n)
        for(j in c:d) {
          result[i, j] <- rss13 * exp(-dists[i, j - b] / rho13)
          result[j, i] <- result[i, j]
        }
      for(i in a:b)
        for(j in c:d) {
          result[i, j] <- rss23 * exp(-dists[i - n, j - b] / rho23)
          result[j, i] <- result[i, j]
        } 
      for(i in a:b)
        for(j in a:b)
          result[i, j] <- sigma2_sq * exp(-dists[i - n, j - n] / rho2)
      for(i in c:d)
        for(j in c:d)
          result[i, j] <- sigma3_sq * exp(-dists[i - b, j - b] / rho3)
      return(result)
    })
  ccovCombine <- compileNimble(covCombine)
  
  # calculate rho_12
  calcRho <- nimbleFunction(
    run = function(rho1 = double(0), rho2 = double(0)) {
      returnType(double(0))
      result = sqrt(2 / ((1 / rho1^2) + (1 / rho2^2)))
      return(result)
    })
  ccalcRho <- compileNimble(calcRho)
  
  # for a user-defined distribution
  assign('expcov', expcov, envir = .GlobalEnv)
  assign('covCombine', covCombine, envir = .GlobalEnv)
  assign('calcRho', calcRho, envir = .GlobalEnv)
  
  # nimble model code
  code <- nimbleCode({
    # global parameters
    log_mean_mu ~ dnorm(0.7, sd = 0.3)
    mu0 ~ dnorm(73, sd = 15)
    mu1 ~ dnorm(-0.3, sd = 0.3)
    
    # MGP parameters
    sigma.a1 ~ T(dnorm(0, sd = 10), 0, Inf)
    sigma.a2 ~ T(dnorm(0, sd = 2.5), 0, Inf)
    sigma.a3 ~ T(dnorm(0, sd = 2.5), 0, Inf)
    sigma.b1 ~ T(dnorm(0, sd = 10), 0, Inf)
    sigma.b2 ~ T(dnorm(0, sd = 2.5), 0, Inf)
    sigma.b3 ~ T(dnorm(0, sd = 2.5), 0, Inf)
    
    rho.a1 ~ dinvgamma(2.25, 2.5)
    rho.a2 ~ dinvgamma(2.25, 2.5)
    rho.a3 ~ dinvgamma(2.25, 2.5)
    rho.a12 <- calcRho(rho1 = rho.a1, rho2 = rho.a2)
    rho.a13 <- calcRho(rho1 = rho.a1, rho2 = rho.a3)
    rho.a23 <- calcRho(rho1 = rho.a2, rho2 = rho.a3)
    S_a[1:3, 1:3] ~ dinvwish(I3[1:3, 1:3], 4)
    R_a12 <- S_a[1, 2] / sqrt(S_a[1, 1] * S_a[2, 2])
    R_a13 <- S_a[1, 3] / sqrt(S_a[1, 1] * S_a[3, 3])
    R_a23 <- S_a[2, 3] / sqrt(S_a[2, 2] * S_a[3, 3])
    
    rho.b1 ~ dinvgamma(2.25, 2.5)
    rho.b2 ~ dinvgamma(2.25, 2.5)
    rho.b3 ~ dinvgamma(2.25, 2.5)
    rho.b12 <- calcRho(rho1 = rho.b1, rho2 = rho.b2)
    rho.b13 <- calcRho(rho1 = rho.b1, rho2 = rho.b3)
    rho.b23 <- calcRho(rho1 = rho.b2, rho2 = rho.b3)
    S_b[1:3, 1:3] ~ dinvwish(I3[1:3, 1:3], 4)
    R_b12 <- S_b[1, 2] / sqrt(S_b[1, 1] * S_b[2, 2])
    R_b13 <- S_b[1, 3] / sqrt(S_b[1, 1] * S_b[3, 3])
    R_b23 <- S_b[2, 3] / sqrt(S_b[2, 2] * S_b[3, 3])
    
    # Visit effects GP parameters
    sigma.v ~ T(dnorm(0, sd = 2.5), 0, Inf)
    rho.v ~ dinvgamma(2.25, 2.5)
    
    # population/superpixel level parameters
    mu[1:36] <- mu0 * ones[1:n_sp]
    mu[37:72] <- mu1 * ones[1:n_sp]
    mu[73:108] <- log_mean_mu * ones[1:n_sp]
    a[1:108] ~ dmnorm(mu[1:108], cov = cov.a[1:108, 1:108])
    
    # covariance matrix for visit effect GP
    cov.v[1:36, 1:36] <- expcov(dists = dists[1:n_sp, 1:n_sp],
                                rho = rho.v, sigma = sigma.v)
    
    # covariance matrices for MGPs
    cov.a[1:108, 1:108] <- covCombine(dists = dists[1:n_sp, 1:n_sp], 
                                      rho1 = rho.a1, rho2 = rho.a2, 
                                      rho3 = rho.a3, rho12 = rho.a12, 
                                      rho13 = rho.a13, rho23 = rho.a23, 
                                      sigma1 = sigma.a1, sigma2 = sigma.a2, 
                                      sigma3 = sigma.a3, R_v12 = R_a12, 
                                      R_v13 = R_a13, R_v23 = R_a23)
    
    cov.b[1:108, 1:108] <- covCombine(dists = dists[1:n_sp, 1:n_sp], 
                                      rho1 = rho.b1, rho2 = rho.b2, 
                                      rho3 = rho.b3, rho12 = rho.b12, 
                                      rho13 = rho.b13, rho23 = rho.b23, 
                                      sigma1 = sigma.b1, sigma2 = sigma.b2, 
                                      sigma3 = sigma.b3, R_v12 = R_b12, 
                                      R_v13 = R_b13, R_v23 = R_b23)
    
    # likelihood
    for(i in 1:n) {
      y[i] ~ dnorm(vis_eff[visit_num[i], sp[i]] + b[id[i], sp[i]] + 
                     (b[id[i], sp2[i]]) * fuyrs[i], sd = exp(b[id[i], sp3[i]]))
    }
    # subject-specific intercepts, slopes, log residual SDs 
    for(j in 1:n_id) {
      b[j, 1:108] ~ dmnorm(a[1:108], cov = cov.b[1:108, 1:108])
    }
    # visit effects
    for(k in 1:n_visits) {
      vis_eff[k, 1:36] ~ dmnorm(zeroes[1:36], cov = cov.v[1:36, 1:36])
    }
  })
  # set up constants and data
  constants <- list(n = n, dists = dists, ones = rep(1, n_sp), n_id = n_id,
                    zeroes = rep(0, n_sp), n_visits = n_visits, 
                    visit_num = visit_num, fuyrs = fuyrs, id = id, sp = sp,
                    n_sp = n_sp, sp2 = sp2, sp3 = sp3, I3 = diag(3))
  data <- list(y = y)
  # initial values
  inits <- list(mu0 = init_data$mu0, 
                mu1 = init_data$mu1, 
                sigma.a1 = init_data$sigma.a1,
                sigma.a2 = init_data$sigma.a2,
                sigma.a3 = init_data$sigma.a3,
                rho.a1 = init_data$rho.a1,
                rho.a2 = init_data$rho.a2,
                rho.a3 = init_data$rho.a3,
                rho.b1 = init_data$rho.b1,
                rho.b2 = init_data$rho.b2,
                rho.b3 = init_data$rho.b3,
                sigma.b1 = init_data$sigma.b1,
                sigma.b2 = init_data$sigma.b2,
                sigma.b3 = init_data$sigma.b3,
                log_mean_mu = init_data$log_mean_mu,
                rho.v = init_data$rho.v,
                sigma.v = init_data$sigma.v,
                S_a = init_data$S_a,
                S_b = init_data$S_b,
                cov.a = init_data$cov.a,
                cov.b = init_data$cov.b,
                cov.v = init_data$cov.v,
                vis_eff = init_data$vis_eff,
                a = init_data$a,
                b = init_data$b
  )
  
  # compile and configure MCMC
  model <- nimbleModel(code, constants = constants, data = data, inits = inits)
  cModel <- compileNimble(model)
  conf <- configureMCMC(model, enableWAIC = T)
  
  # monitor all parameters of interest
  conf$addMonitors(c("a", "b", "rho.a12", "rho.a13", "rho.a23", "rho.b12", 
                     "rho.b13", "rho.b23", "R_a12", "R_a13", "R_a23", "R_b12", 
                     "R_b13", "R_b23", "vis_eff"))
  # customize samplers
  conf$removeSampler(c("S_a", "S_b"))
  conf$addSampler("S_a[1:3, 1:3]", 'RW_wishart',
                  control = list(adaptFactorExponent = 0.8, 
                                 adaptInterval = 200), silent = T)   
  conf$addSampler("S_b[1:3, 1:3]", 'RW_wishart',
                  control = list(adaptFactorExponent = 0.8, 
                                 adaptInterval = 200), silent = T)   
  conf$removeSampler(c("b[]"))
  conf$removeSampler(c("vis_eff[]"))
  
  first_visits = n_i$visit_num[n_i$visit == 1]
  sampler_nums = seq(0, 36, 3)
  for (i in 1:n_id) {
    for (j in 1:12) {
      conf$addSampler(c(paste0("b[", i, ", ", sampler_nums[j] + 1, ":",
                               sampler_nums[j + 1], "]"),
                        paste0("b[", i, ", ", sampler_nums[j] + 37, ":",
                               sampler_nums[j + 1] + 36, "]"),
                        paste0("vis_eff[", first_visits[i], ", ", 
                               sampler_nums[j] + 1, ":",
                               sampler_nums[j + 1], "]")), 'RW_block',
                      control = list(adaptFactorExponent = 0.8, 
                                     adaptInterval = 200), silent = T)   
    }
  }
  
  sampler_nums = seq(72, 108, 6)
  for (i in 1:n_id) {
    for (j in 1:6) {
      conf$addSampler(paste0("b[", i, ", ", sampler_nums[j] + 1, ":",
                             sampler_nums[j + 1], "]"), 'RW_block',
                      control = list(adaptFactorExponent = 0.8, 
                                     adaptInterval = 200), silent = T)   
    }
  }
  
  visit_list = seq(1, n_visit)
  visit_list = visit_list[!visit_list %in% first_visits]
  sampler_nums = seq(0, 36, 3)
  for (i in 1:length(visit_list)) {
    for (j in 1:12) {
      conf$addSampler(paste0("vis_eff[", visit_list[i], ", ", 
                             sampler_nums[j] + 1, ":",
                             sampler_nums[j + 1], "]"), 'RW_block',
                      control = list(adaptFactorExponent = 0.8, 
                                     adaptInterval = 20), silent = T)   
    }
  }
  
  conf$removeSamplers('sigma.b1', 'rho.b1', 'sigma.b2', 'rho.b2', 'sigma.b3', 
                      'rho.b3')
  conf$addSampler(c('sigma.b1', 'rho.b1'), 'RW_block',
                  control = list(adaptFactorExponent = 0.8, adaptInterval = 20),
                  silent = T)
  conf$addSampler(c('sigma.b2', 'rho.b2'), 'RW_block',
                  control = list(adaptFactorExponent = 0.8, adaptInterval = 20),
                  silent = T)
  conf$removeSamplers('sigma.a1', 'rho.a1', 'sigma.a2', 'rho.a2', 'sigma.a3', 
                      'rho.a3')
  conf$addSampler(c('sigma.a1', 'rho.a1'), 'RW_block',
                  control = list(adaptFactorExponent = 0.8, adaptInterval = 20),
                  silent = T)
  conf$addSampler(c('sigma.a2', 'rho.a2'), 'RW_block',
                  control = list(adaptFactorExponent = 0.8, adaptInterval = 20),
                  silent = T)
  conf$addSampler(c('sigma.a3', 'rho.a3'), 'RW_block',
                  control = list(adaptFactorExponent = 0.8, adaptInterval = 20),
                  silent = T)
  conf$removeSamplers('sigma.v', 'rho.v')
  conf$addSampler(c('sigma.v', 'rho.b3'), 'RW_block',
                  control = list(adaptFactorExponent = 0.8, adaptInterval = 20),
                  silent = T)
  conf$addSampler(c('rho.v', 'sigma.b3'), 'AF_slice',
                  control = list(sliceAdaptFactorInterval = 1000),
                  silent = T)
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC, project = cModel)
  # specify iterations, burn-in, and thinning for each chain
  results <- runMCMC(cMCMC, niter = 250000, setSeed = seed, nburnin = 40000,
                     thin = 100)
  return(results)
}

# ---- Parallel Run ----
# running chains in parallel (X is number of chains)
set.seed(502)
chain_output <- parLapply(cl = this_cluster, X = 1:9, 
                          fun = run_MCMC_allcode, 
                          SP_data = SP_data, init_data = init_data)

# Close cluster
stopCluster(this_cluster)

# save model file as RDS
saveRDS(chain_output, file = "model.rds")
