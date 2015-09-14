#!/usr/bin/Rscript --vanilla
#
# evalute_results
# %DESCRIPTION%
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-09-05
# Licence:  GPLv3 or later
#

res.version <- c()
res.model   <- c()
res.mse     <- c()
res.run.time <- c()

options(digits=2)

versions <- list.dirs("runs", recursive=F)
for (version in versions) {
  models <- list.dirs(version, recursive=F)
  for (model in models) {
    file.est <- paste(model, 'results', "estimates.txt", sep="/")
    file.true.values <- paste(model, 'results', "true_values.txt", sep="/")
    file.run.times <- paste(model, 'results', "runtimes.txt", sep="/")

    if (!(file.exists(file.est) & file.exists(file.true.values) &
          file.exists(file.run.times)) ) {
      warning(version, " ", model, ": Some results are missing")
      next()
    }

    estimates <- read.table(file.est, header=T)
    true.values <- read.table(file.true.values, header=T)
    run.times <- read.table(file.run.times, header=T)
    if (any(is.na(estimates))) {
      warning(version, " ", model, ": Some results are missing")
      next()
    }

    par.error <- apply(abs(log10(estimates/true.values)), 2, mean)
    print(model)
    print(par.error)
    cat('\n') 
    mse <- mean(par.error)
    if (ncol(run.times) == 6) {
      avg_run_time <- mean(run.times[, 3]) +  mean(run.times[, 6])
    } else {
      avg_run_time <- mean(run.times[, 3])
    }

    version <- strsplit(model, "/")[[1]][2]
    model <- strsplit(model, "/")[[1]][3]

    res.version <- c(res.version, version)
    res.model   <- c(res.model, model)
    res.mse     <- c(res.mse, mse)
    res.run.time <- c(res.run.time, avg_run_time)
  }
}

result <- data.frame(Version = res.version, 
           Model = res.model,
           Error = res.mse,
           Run.Time = res.run.time)

result <- result[with(result, order(Model, Version)), ]

print(result)
