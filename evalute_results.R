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

versions <- list.dirs("results", recursive=F)
for (version in versions) {
  version.name <- list.dirs()
  models <- list.dirs(version, recursive=F)
  for (model in models) {
    estimates <- read.table(paste(model, "/estimates.txt", sep=""), header=T)
    true.values <- read.table(paste(model, "/true_values.txt", sep=""), header=T)
    run.times <- read.table(paste(model, "/runtimes.txt", sep=""), header=T)

    mse <- mean((estimates - true.values)^2)
    avg.run.time <- mean(run.times[, 3] + run.times[,6])
    version <- strsplit(models, "/")[[1]][2]
    model <- strsplit(models, "/")[[1]][3]

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

result
