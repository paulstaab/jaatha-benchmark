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

library(reshape2)
library(foreach)

list.dirs('runs')

runs <- lapply(list.dirs("runs", recursive=FALSE), 
               function(x) list.dirs(x, recursive=FALSE))

data_acc <- foreach(model=unlist(runs), .combine=rbind) %do% {
  file.est <- paste(model, 'results', "estimates.txt", sep="/")
  file.true.values <- paste(model, 'results', "true_values.txt", sep="/")
  file.run.times <- paste(model, 'results', "runtimes.txt", sep="/")
  
  if (!(file.exists(file.est) & file.exists(file.true.values) &
          file.exists(file.run.times)) ) {
    warning(version, " ", model, ": Some results are missing")
    return()
  }
  
  estimates <- read.table(file.est, header=T)
  true.values <- read.table(file.true.values, header=T)
  run.times <- read.table(file.run.times, header=T)
  
  version_name <- strsplit(model, "/")[[1]][2]
  model_name <- strsplit(model, "/")[[1]][3]
  
  rel.error <- (true.values-estimates)/true.values
  colnames(rel.error) <- colnames(estimates)
  cbind(version=version_name, 
        model=model_name,
        melt(rel.error, variable.names="parameter", value.name="rel.error"))
}

data_runtime <- foreach(model=unlist(runs), .combine=rbind) %do% {
  file.est <- paste(model, 'results', "estimates.txt", sep="/")
  file.true.values <- paste(model, 'results', "true_values.txt", sep="/")
  file.run.times <- paste(model, 'results', "runtimes.txt", sep="/")
  
  if (!(file.exists(file.est) & file.exists(file.true.values) &
          file.exists(file.run.times)) ) {
    warning(version, " ", model, ": Some results are missing")
    return()
  }
  
  estimates <- read.table(file.est, header=T)
  true.values <- read.table(file.true.values, header=T)
  run.times <- read.table(file.run.times, header=T)
  
  version_name <- strsplit(model, "/")[[1]][2]
  model_name <- strsplit(model, "/")[[1]][3]
  
  runtimes <- run.times[ , 3]
  if (ncol(run.times) == 10) {
    runtimes <- runtimes + run.times[ , 8]
  }
    
  data.frame(version=version_name,
             model=model_name,
             run.time=runtimes)
}

save(data_acc, data_runtime, file='data_processed/acc_and_runtime.Rda')
