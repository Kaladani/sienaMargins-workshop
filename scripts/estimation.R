##############################################
# R script analyzing AK friendship data
# for waves W -> X in all available classes;
# Script 1: Estimation -- all six specifications.
# version of 02 July 2026 (models 0-2) and 17 July 2026 (models 3-5)
# written by Daniel Gotthardt, Christian Steglich & Marijtie van Duijn
##############################################
#
# Written to be read and stepped through: one block per specification, each
# followed by its convergence check and, where classes needed a second pass,
# the rerun that fixed them.  The notes about which classes those were are
# left in place -- they are what running this actually looks like.
#
# Everything is written once at the end to
#
#     results/estimation_results.rds
#
# which is the single file every later script reads.
#
# To adapt this to your own data, change three things: the data block below,
# the model_N specification blocks, and the `selection` vector naming your
# groups.  Each model is estimated in its own block, so you can run one and
# stop.
#
# This is the slow step -- hours for all six specifications across 17 classes.
# The output ships with this bundle, so run it only if you change the data or
# the specifications.  For a smaller worked version end to end, see
# 04-minimal-meta-analysis.qmd.


#Dependency

# using marginalEffects branch
# renv::install("stocnet/rsiena@sienaMargins")
# OR
# devtools::install_github("stocnet/rsiena@sienaMargins")

library(RSiena)
library(sna)

# DATA IMPORT AND FIRST INSPECTION ----
## Change if necessary
data_dir <- "data/"

## The AK data are not distributed with this material.  This fetches them into
## data_dir on first use; if they are already there it does nothing.  See
## data/README.md for the source and the terms.
source("scripts/download_data.R")
ensure_ak_data(dir = data_dir)

load(file.path(data_dir, "AK_friendship.RData"))
load(file.path(data_dir, "AK_gender.RData"))
load(file.path(data_dir, "AK_primary.RData"))

classIDs <- names(friendship)
names(classIDs) <- classIDs

school12_names <- c("12b", "12c", "12e", "12f", "12g", "12h", "12k", "12m", "12n", "12p")
school02_names <- c("02a", "02b", "02c", "02d", "02e", "02f", "02g")
selection <- c(school12_names, school02_names)


# RSIENA ANALYSES (NEW SYNTAX) -----
thedata <- lapply(selection, function(cl){
  numberActors <- length(male[[cl]])
  dat <- make_data_rsiena(
    friendship = as_dependent_rsiena(array(
      c(friendship[[cl]][['W']], friendship[[cl]][['X']]),
      dim=c(numberActors, numberActors, 2)
    )),
    sex_m = as_covariate_rsiena(male[[cl]], centered=FALSE),
    primary = as_covariate_rsiena(primary[[cl]], type="oneMode", centered=FALSE)
  )
  return(dat)
})


# ---- specifications -----------------------------------------------------
model_0 <- lapply(thedata, function(td){
  model <- make_specification(td)
  model <- set_effect(model, c(recip, transTrip1), verbose=FALSE)
  model <- set_interaction(model, c(recip, transTrip1))
  return(model)
})



model_1 <- lapply(thedata, function(td){
  model <- make_specification(td)
  model <- set_effect(model, c(egoX, altX, egoXaltX), covar1='sex_m', verbose=FALSE)
  model <- set_effect(model, X, covar1='primary', verbose=FALSE)
  model <- set_effect(model, recip, include=FALSE, verbose=FALSE)
  return(model)
})



model_2 <- lapply(model_1, function(m1){
  model <- set_effect(m1, c(recip, transTrip1), verbose=FALSE)
  model <- set_interaction(model, c(recip, transTrip1))
  return(model)
})


model_3 <- lapply(model_1, function(m1){
  model <- set_effect(m1, c(recip, transTrip1), verbose=FALSE)
  return(model)
})


model_4 <- lapply(thedata, function(td){
  model <- make_specification(td)
  model <- set_effect(model, c(egoX, altX, egoXaltX), covar1='sex_m', verbose=FALSE)
  model <- set_effect(model, X, covar1='primary', verbose=FALSE)
  return(model)
})


model_5 <- lapply(thedata, function(td){
  model <- make_specification(td)
  model <- set_effect(model, c(egoX, altX), covar1='sex_m', verbose=FALSE)
  model <- set_effect(model, X, covar1='primary', verbose=FALSE)
  model <- set_effect(model, recip, include=TRUE, verbose=FALSE)
  return(model)
})


# ---- algorithm settings -------------------------------------------------
alg_controls <- lapply(selection, function(cl){
  controls <- set_algorithm_saom(cond = FALSE, n3 = 1500, seed = 1234567)
  return(controls)
})


alg_controls_rerun <- lapply(selection, function(cl){
  controls <- set_algorithm_saom(cond = FALSE, nsub = 1, n2start = 250, n3 = 4000, firstg = 0.02, seed = 1234567)
  return(controls)
})

out_controls <- lapply(selection, function(cl){
  controls <- set_output_saom(outputName = paste0('class_', cl))
  return(controls)
})
# ---- estimation ---------------------------------------------------------
# estimate the structural model:
results_0 <- list()
for (i in 1:length(selection)){
  cat('\r', paste0('Processing class_', selection[i], ' (', i,'/', length(selection), ').\n'))
  flush.console()
  results_0[[i]] <- siena(data = thedata[[i]], effects = model_0[[i]], 
                          control_algo = alg_controls[[i]], control_out = out_controls[[i]], 
                          batch = TRUE, verbose = FALSE)
}
names(results_0) <- selection

# check convergence indicators:
(problem_0 <- which(sapply(results_0, function(rs){
  return(rs$tconv.max >= .25 | rs$tmax >= .1)
})))

# 12e 12f are problematic

for (i in problem_0){
  cat('\r', paste0('Rerunning class_', selection[i], ' (', which(problem_0==i),'/', length(problem_0), ').\n'))
  flush.console()
  results_0[[i]] <- siena(data = thedata[[i]], effects = model_0[[i]], 
                          control_algo = alg_controls_rerun[[i]], control_out = out_controls[[i]], 
                          batch = TRUE, verbose = FALSE, prevAns = results_0[[i]])
}

# check convergence again:
(problem_0 <- which(sapply(results_0, function(rs){
  return(rs$tconv.max >= .25 | rs$tmax >= .1)
})))

# no issues left



# estimate the first model:
results_1 <- list()
for (i in 1:length(selection)){
  cat('\r', paste0('Processing class_', selection[i], ' (', i,'/', length(selection), ').\n'))
  flush.console()
  results_1[[i]] <- siena(data = thedata[[i]], effects = model_1[[i]], 
                          control_algo = alg_controls[[i]], control_out = out_controls[[i]], 
                          batch = TRUE, verbose = FALSE)
}
names(results_1) <- selection

# check convergence indicators:
(problem_1 <- which(sapply(results_1, function(rs){
  return(rs$tconv.max >= .25 | rs$tmax >= .1)
})))


# estimate the second model:
results_2 <- list()
for (i in 1:length(selection)){
  cat('\r', paste0('Processing class_', selection[i], ' (', i,'/', length(selection), ').\n'))
  flush.console()
  results_2[[i]] <- siena(data = thedata[[i]], effects = model_2[[i]], 
                          control_algo = alg_controls[[i]], control_out = out_controls[[i]], 
                          batch = TRUE, verbose = FALSE)
}
names(results_2) <- selection

# check convergence and rerun
(problem_2 <- which(sapply(results_2, function(rs){ return(rs$tconv.max >= .25 | rs$tmax >= .1) 
})))
# 12 c and 02e have issues but only barely
for (i in problem_2){
  cat('\r', paste0('Rerunning class_', selection[i], ' (', which(problem_2==i),'/', length(problem_2), ').\n'))
  flush.console()
  results_2[[i]] <- siena(data = thedata[[i]], effects = model_2[[i]], 
                          control_algo = alg_controls_rerun[[i]], control_out = out_controls[[i]], 
                          batch = TRUE, verbose = FALSE, prevAns = results_2[[i]])
}

# check convergence indicators:
(problem_2 <- which(sapply(results_2, function(rs){
  return(rs$tconv.max >= .25 | rs$tmax >= .1)
})))
# no problems left

# estimate model 3:
results_3 <- list()
for (i in 1:length(selection)){
  cat('\r', paste0('Processing class_', selection[i], ' (', i,'/', length(selection), ').\n'))
  flush.console()
  results_3[[i]] <- siena(data = thedata[[i]], effects = model_3[[i]], 
                          control_algo = alg_controls[[i]], control_out = out_controls[[i]], 
                          batch = TRUE, verbose = FALSE)
}
names(results_3) <- selection

# check convergence indicators:
(problem_3 <- which(sapply(results_3, function(rs){
  return(rs$tconv.max >= .25 | rs$tmax >= .1)
})))

# no issues


# estimate model 4
results_4 <- list()
for (i in 1:length(selection)){
  cat('\r', paste0('Processing class_', selection[i], ' (', i,'/', length(selection), ').\n'))
  flush.console()
  results_4[[i]] <- siena(data = thedata[[i]], effects = model_4[[i]], 
                          control_algo = alg_controls[[i]], control_out = out_controls[[i]], 
                          batch = TRUE, verbose = FALSE)
}
names(results_4) <- selection

# check convergence indicators:
(problem_4 <- which(sapply(results_4, function(rs){
  return(rs$tconv.max >= .25 | rs$tmax >= .1)
})))

# no issues


# estimate model 5:
results_5 <- list()
for (i in 1:length(selection)){
  cat('\r', paste0('Processing class_', selection[i], ' (', i,'/', length(selection), ').\n'))
  flush.console()
  results_5[[i]] <- siena(data = thedata[[i]], effects = model_5[[i]], 
                          control_algo = alg_controls[[i]], control_out = out_controls[[i]], 
                          batch = TRUE, verbose = FALSE)
}
names(results_5) <- selection

# check convergence indicators:
(problem_5 <- which(sapply(results_5, function(rs){
  return(rs$tconv.max >= .25 | rs$tmax >= .1)
})))

# no issues

# ---- one combined save ---------------------------------------------------

# Compile Raw Coefficients
extract_raw_coefs <- function(res_list, model_name) {
  do.call(rbind, lapply(names(res_list), function(cl) {
    fit <- res_list[[cl]]
    
    # Force RSiena to return the rate parameters alongside the structural effects
    estimates <- coef(fit, dropRates = FALSE)
    standard_errors <- sqrt(diag(vcov(fit, dropRates = FALSE)))
    
    data.frame(
      est      = estimates, 
      se       = standard_errors,
      effect   = names(estimates), 
      class_id = cl, 
      model    = model_name, 
      row.names = NULL
    )
  }))
}

# Each model paired with its own fits.  Worth checking by eye when you add a
# specification: a mismatch here would relabel one model's coefficients as
# another's, and nothing downstream could detect it.
results <- list(model_0 = results_0, model_1 = results_1, model_2 = results_2,
                model_3 = results_3, model_4 = results_4, model_5 = results_5)

models <- list(model_0 = model_0, model_1 = model_1, model_2 = model_2,
               model_3 = model_3, model_4 = model_4, model_5 = model_5)

problems <- list(model_0 = problem_0, model_1 = problem_1, model_2 = problem_2,
                 model_3 = problem_3, model_4 = problem_4, model_5 = problem_5)

coefs <- do.call(rbind, lapply(names(results), function(nm)
  extract_raw_coefs(results[[nm]], nm)))

convergence <- data.frame(
  model         = names(results),
  not_converged = vapply(problems, function(p)
                    paste(names(p), collapse = " "), character(1L)),
  row.names = NULL)

names(thedata)      <- selection
names(alg_controls) <- selection

estimation <- list(selection    = selection,
                   thedata      = thedata,
                   models       = models,
                   alg_controls = alg_controls,
                   results      = results,
                   coefs        = coefs,
                   convergence  = convergence,
                   when         = Sys.time())

saveRDS(estimation, "results/estimation_results.rds")

