#!/usr/bin/Rscript --vanilla
#
# evalute_results
# %DESCRIPTION%
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-27
# Licence:  GPLv3 or later
#

args <- commandArgs(TRUE)

# Enable just-in-time-compiler if availible
if ("compiler" %in% rownames(installed.packages())){
  library("compiler")
  invisible(compiler::enableJIT(3))
}

load("./testResults.save")
if (!exists("test.results")) stop("No results found!")

calc.mse <- function(test.results) {
  res.version <- c()
  res.model   <- c()
  res.mse     <- c()
  res.run.time <- c()

  for (i in seq(along = test.results)){
    version <- names(test.results)[i]
    for (j in seq(along = test.results[[i]])) {
      model <- names(test.results[[i]])[j]
      current.result <- test.results[[i]][[j]]
      mse <- mean((current.result$RealPars - current.result$Estimates)^2)
      avg.run.time <- mean(current.result$RunTimes[ ,8])
      res.version <- c(res.version, version)
      res.model   <- c(res.model, model)
      res.mse     <- c(res.mse, mse)
      res.run.time <- c(res.run.time, avg.run.time)
    }
  }

  return(data.frame(version=res.version, model=res.model, 
                    mse=res.mse, avg.run.time=res.run.time))
}

calc.mse(test.results)
