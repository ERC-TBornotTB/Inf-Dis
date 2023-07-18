# Functions to run cohort simulation

full_simulation_inf <- function(n_steps, n_people, prop_inf, params, method, treatment=0, prop_min=0, prop_sub=0){
  pop <- initialise_population_inf(n_steps, n_people, prop_inf, prop_min, prop_sub)
  risks <- sort_risks_inf(n_steps, n_people, params, method)
  for (i in seq(2,n_steps)) {
    pop <- update_population_inf(pop, i, n_people, risks, treatment)
  }
  return(pop)
}

initialise_population_inf <- function(n_steps, n_people, prop_inf, prop_min=0, prop_sub=0){
  if (prop_inf < 0 || prop_inf > 1) stop("prop_inf is out of the range 0 - 1")
  if (prop_inf + prop_sub + prop_min > 1) stop("prop people over 1")
  n_inf <- round(prop_inf * n_people)
  n_min <- round(prop_min * n_people)
  n_sub <- round(prop_sub * n_people)
  n_clin <- n_people - (n_sub + n_min + n_inf)
  I <- matrix(data = NA, nrow = n_steps, ncol = n_people)
  I[1,] <- c(rep(x = "i", times = n_inf), rep(x = "m", times = n_min), rep(x = "s", times = n_sub), rep(x = "c", times = n_clin))
  return(I)
}

param_to_risk <- function(param){
  param_new <- 1 - exp(-param/12)
  return(param_new)
}

sort_risks_inf <- function(n_steps, n_people, params, method){
  if (isFALSE(dplyr::setequal(names(params), c("inf_clr","inf_sub","inf_min","min_rec","min_sub","sub_min","sub_clin","clin_sub","mort")))) stop("incorrect parameter names - include all of min_sub, sub_min, sub_clin, clin_sub")
  if (isFALSE(method %in% c("median","fixed","random"))) stop("incorrect method - choose one of median, fixed, or random")
  params <- as.data.frame(params)
  if (method == "median") {
    R = median_risks_inf(n_steps, n_people, params)
  } else if (method == "fixed") {
    R = fixed_risks_inf(n_steps, n_people, params)
  } else if (method == "random") {
    R = random_risks_inf(n_steps, n_people, params)
  }
  return(R)
}        

median_risks_inf <- function(n_steps, n_people, params){
  med_params <- apply(params, MARGIN = 2, FUN = median)
  med_risks <- param_to_risk(med_params)
  R <- list(inf_clr  = matrix(data = med_risks[["inf_clr"]],  nrow = n_steps, ncol = n_people),
            inf_sub  = matrix(data = med_risks[["inf_sub"]],  nrow = n_steps, ncol = n_people),
            inf_min  = matrix(data = med_risks[["inf_min"]],  nrow = n_steps, ncol = n_people),
            min_rec  = matrix(data = med_risks[["min_rec"]],  nrow = n_steps, ncol = n_people),
            min_sub  = matrix(data = med_risks[["min_sub"]],  nrow = n_steps, ncol = n_people),
            sub_min  = matrix(data = med_risks[["sub_min"]],  nrow = n_steps, ncol = n_people),
            sub_clin = matrix(data = med_risks[["sub_clin"]], nrow = n_steps, ncol = n_people),
            clin_sub = matrix(data = med_risks[["clin_sub"]], nrow = n_steps, ncol = n_people),
            mort     = matrix(data = med_risks[["mort"]],     nrow = n_steps, ncol = n_people))
  return(R)
}

fixed_risks_inf <- function(n_steps, n_people, params){
  risks <- param_to_risk(params)
  R <- list(inf_clr  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            inf_sub  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            inf_min  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            min_rec  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            min_sub  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            sub_min  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            sub_clin = matrix(data = NA, nrow = n_steps, ncol = n_people),
            clin_sub = matrix(data = NA, nrow = n_steps, ncol = n_people),
            mort     = matrix(data = NA, nrow = n_steps, ncol = n_people))
  for (i in 1:n_people) {
    R$inf_clr[,i]  = rep(base::sample(x = risks$inf_clr,  size = 1), n_steps)
    R$inf_sub[,i]  = rep(base::sample(x = risks$inf_sub,  size = 1), n_steps)
    R$inf_min[,i]  = rep(base::sample(x = risks$inf_min,  size = 1), n_steps)
    R$min_rec[,i]  = rep(base::sample(x = risks$min_rec,  size = 1), n_steps)
    R$min_sub[,i]  = rep(base::sample(x = risks$min_sub,  size = 1), n_steps)
    R$sub_min[,i]  = rep(base::sample(x = risks$sub_min,  size = 1), n_steps)
    R$sub_clin[,i] = rep(base::sample(x = risks$sub_clin, size = 1), n_steps)
    R$clin_sub[,i] = rep(base::sample(x = risks$clin_sub, size = 1), n_steps)
    R$mort[,i]     = rep(base::sample(x = risks$mort,     size = 1), n_steps)
  }
  return(R)
}

random_risks_inf <- function(n_steps, n_people, params){
  risks <- param_to_risk(params)
  R <- list(inf_clr  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            inf_sub  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            inf_min  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            min_rec  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            min_sub  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            sub_min  = matrix(data = NA, nrow = n_steps, ncol = n_people),
            sub_clin = matrix(data = NA, nrow = n_steps, ncol = n_people),
            clin_sub = matrix(data = NA, nrow = n_steps, ncol = n_people),
            mort     = matrix(data = NA, nrow = n_steps, ncol = n_people))
  for (i in 1:n_people) {
    R$inf_clr[,i]  = base::sample(x = risks$inf_clr,  size = n_steps, replace = TRUE)
    R$inf_sub[,i]  = base::sample(x = risks$inf_sub,  size = n_steps, replace = TRUE)
    R$inf_min[,i]  = base::sample(x = risks$inf_min,  size = n_steps, replace = TRUE)
    R$min_rec[,i]  = base::sample(x = risks$min_rec,  size = n_steps, replace = TRUE)
    R$min_sub[,i]  = base::sample(x = risks$min_sub,  size = n_steps, replace = TRUE)
    R$sub_min[,i]  = base::sample(x = risks$sub_min,  size = n_steps, replace = TRUE)
    R$sub_clin[,i] = base::sample(x = risks$sub_clin, size = n_steps, replace = TRUE)
    R$clin_sub[,i] = base::sample(x = risks$clin_sub, size = n_steps, replace = TRUE)
    R$mort[,i]     = base::sample(x = risks$mort,     size = n_steps, replace = TRUE)
  }
  return(R)
}

update_population_inf <- function(X, t, n_people, risks, treatment){
  for (i in 1:n_people) {
    if (X[t - 1,i] == "i") {
      X[t,i] <-  update_infect_inf(i, t-1, risks)
    } else if (X[t - 1,i] == "m") {
      X[t,i] <-  update_minimal_inf(i, t-1, risks)
    } else if (X[t - 1,i] == "s") {
      X[t,i] <-  update_subclinical_inf(i, t-1, risks)
    } else if (X[t - 1,i] == "c") {
      X[t,i] <-  update_clinical_inf(i, t-1, risks, treatment)
    } else {
      X[t,i] <-  "-"
    }
  }
  return(X)
}

update_infect_inf <- function(person, time, risks, treatment){
  inf_sub <- risks$inf_sub[time, person]
  inf_min <- risks$inf_min[time, person]
  inf_clr <- risks$inf_clr[time, person]
  prob <- runif(1,0,1)
  if (prob < inf_clr) {
    return("x")
  } else if (prob < inf_clr + inf_min) {
    return("m")
  } else if (prob < inf_clr + inf_min + inf_sub) {
    return("s")
  } else {
    return("i")
  }
}

update_minimal_inf <- function(person, time, risks){
  min_rec <- risks$min_rec[time, person]
  min_sub <- risks$min_sub[time, person]
  prob <- runif(1,0,1)
  if (prob < min_rec) {
    return("r")
  } else if (prob < min_rec + min_sub) {
    return("s")
  } else {
    return("m")
  }
}

update_subclinical_inf <- function(person, time, risks){
  sub_min  <- risks$sub_min[time, person]
  sub_clin <- risks$sub_clin[time, person]
  prob <- runif(1,0,1)
  if (prob < sub_min) {
    return("m")
  } else if (prob < sub_min + sub_clin) {
    return("c")
  } else {
    return("s")
  }
}

update_clinical_inf <- function(person, time, risks, treatment){
  treat    <- param_to_risk(treatment)
  mort     <- risks$mort[time, person]
  clin_sub <- risks$clin_sub[time, person]
  prob <- runif(1,0,1)
  if (prob < mort) {
    return("d")
  } else if (prob < mort + clin_sub) {
    return("s")
  } else {
    prob_treat <- runif(1,0,1)
    if (prob_treat < treat) {
      return("t")
    } else {
      return("c") 
    }
  }
}

# Functions to calculate incidence

full_simulation_inc <- function(n_steps, n_people, prop_inf, params, method, endpoint, treatment=0, prop_min=0, prop_sub=0){
  pop <- initialise_population_inf(n_steps, n_people, prop_inf, prop_min, prop_sub)
  risks <- sort_risks_inf(n_steps, n_people, params, method)
  for (i in seq(2,n_steps)) {
    pop <- update_population_inc(pop, i, n_people, risks, treatment, endpoint)
  }
  return(pop)
}

update_population_inc <- function(X, t, n_people, risks, treatment, endpoint){
  for (i in 1:n_people) {
    if (endpoint == "clin"){
      if (X[t - 1,i] == "i") {
        X[t,i] <-  update_infect_incclin(i, t-1, risks)
      } else if (X[t - 1,i] == "m") {
        X[t,i] <-  update_minimal_incclin(i, t-1, risks)
      } else if (X[t - 1,i] == "s") {
        X[t,i] <-  update_subclinical_incclin(i, t-1, risks)
      } else {
        X[t,i] <-  "-"
      }
    } else if (endpoint == "sub"){
      if (X[t - 1,i] == "i") {
        X[t,i] <-  update_infect_incsub(i, t-1, risks)
      } else if (X[t - 1,i] == "m") {
        X[t,i] <-  update_minimal_incsub(i, t-1, risks)
      } else {
        X[t,i] <-  "-"
      }
    } else if (endpoint == "min"){
      if (X[t - 1,i] == "i") {
        X[t,i] <-  update_infect_incmin(i, t-1, risks)
      } else {
        X[t,i] <-  "-"
      }
    }
  }
  return(X)
}

update_infect_incclin <- function(person, time, risks, treatment){
  inf_sub <- risks$inf_sub[time, person]
  inf_min <- risks$inf_min[time, person]
  inf_clr <- risks$inf_clr[time, person]
  prob <- runif(1,0,1)
  if (prob < inf_clr) {
    return("r")
  } else if (prob < inf_min + inf_clr) {
    return("m")
  } else if (prob < inf_sub + inf_min + inf_clr) {
    return("s")
  } else {
    return("i")
  }
}
update_minimal_incclin <- function(person, time, risks){
  min_rec <- risks$min_rec[time, person]
  min_sub <- risks$min_sub[time, person]
  prob <- runif(1,0,1)
  if (prob < min_rec) {
    return("r")
  } else if (prob < min_sub + min_rec) {
    return("s")
  } else {
    return("m")
  }
}
update_subclinical_incclin <- function(person, time, risks){
  sub_min  <- risks$sub_min[time, person]
  sub_clin <- risks$sub_clin[time, person]
  prob <- runif(1,0,1)
  if (prob < sub_min) {
    return("m")
  } else if (prob < sub_clin + sub_min) {
    return("c")
  } else {
    return("s")
  }
}
update_clinical_incclin <- function(person, time, risks, treatment){
  return("c") 
}

update_infect_incsub <- function(person, time, risks, treatment){
  inf_sub <- risks$inf_sub[time, person]
  inf_min <- risks$inf_min[time, person]
  inf_clr <- risks$inf_clr[time, person]
  prob <- runif(1,0,1)
  if (prob < inf_clr) {
    return("r")
  } else if (prob < inf_min + inf_clr) {
    return("m")
  } else if (prob < inf_sub + inf_min + inf_clr) {
    return("s")
  } else {
    return("i")
  }
}
update_minimal_incsub <- function(person, time, risks){
  min_rec <- risks$min_rec[time, person]
  min_sub <- risks$min_sub[time, person]
  prob <- runif(1,0,1)
  if (prob < min_rec) {
    return("r")
  } else if (prob < min_sub + min_rec) {
    return("s")
  } else {
    return("m")
  }
}
update_subclinical_incsub <- function(person, time, risks){
  return("s")
}

update_infect_incmin <- function(person, time, risks, treatment){
  inf_sub <- risks$inf_sub[time, person]
  inf_min <- risks$inf_min[time, person]
  inf_clr <- risks$inf_clr[time, person]
  prob <- runif(1,0,1)
  if (prob < inf_clr) {
    return("r")
  } else if (prob < inf_min + inf_clr) {
    return("m")
  } else if (prob < inf_sub + inf_min + inf_clr) {
    return("s")
  } else {
    return("i")
  }
}
update_minimal_incmin <- function(person, time, risks){
  return("m")
}