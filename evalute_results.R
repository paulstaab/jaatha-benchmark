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

args <- commandArgs(TRUE)

# Enable just-in-time-compiler if availible
if ("compiler" %in% rownames(installed.packages())){
  library("compiler")
  invisible(compiler::enableJIT(3))
}

res.version <- c()
res.model   <- c()
res.mse     <- c()
res.run.time <- c()

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

    #mse <- mean((estimates - true.values)^2)
    mse <- mean(apply(abs(log10(estimates/true.values)), 2, mean))
    avg.run.time <- mean(run.times[, 3] + run.times[,6])
    version <- strsplit(model, "/")[[1]][2]
    model <- strsplit(model, "/")[[1]][3]

    res.version <- c(res.version, version)
    res.model   <- c(res.model, model)
    res.mse     <- c(res.mse, mse)
    res.run.time <- c(res.run.time, avg.run.time)
  }
}

result <- data.frame(Version=res.version, 
           Model=res.model,
           MSE=res.mse,
           Run.Time=res.run.time )

print(result)
