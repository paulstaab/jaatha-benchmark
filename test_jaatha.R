#!/usr/bin/Rscript --vanilla

# Tests the current installed version of Jaatha using
# a simple theta/tau model

library(jaatha)
library(foreach)
library(doMC)

result.folder <- "results"

version <- as.character(packageVersion("jaatha"))

cat("Testing version",version,"of Jaatha\n")

args <- commandArgs(TRUE)
if (is.na(args[1])) args[1] <- 1
threads <- as.numeric(args[1])
cat("Using", threads, "Cores\n")

registerDoMC(threads)

runTest <- function(dm, n.points=5, seed=12523, model){
  set.seed(seed)

  par.grid <- CreateParGrid(dm, n.points)

  cat("Simulating data...\n")
  sumStats <- dm.simSumStats(dm,par.grid)

  n <- dim(par.grid)[1]
  seeds <- sample.int(2^20, n)


  folder <- paste("results", "/", version, "/", model, "/", sep="")
  log.folder <- paste("logs", "/", version, "/", model, "/", sep="")
  dir.create(folder, recursive=T, showWarnings=F)
  dir.create(log.folder, recursive=T, showWarnings=F)

  results <- foreach(i=1:n, .combine=rbind) %dopar% { 
      cat("Run",i,"of",n,"\n")
      sink(file=paste(log.folder, "run_", i, ".txt", sep=""))
      set.seed(seeds[i])
      cat("----------------------------------------------------------------------\n")
      cat("Run",i,"of",n,"\n")
      cat("Real parameters:",par.grid[i,],"\n")
      cat("----------------------------------------------------------------------\n")
      jaatha <- Jaatha.initialize(dm, sumStats[i,])

      runtimes <- rep(0, 6)
      names(runtimes) <-
        c('init.user','init.system','init.elapsed','ref.user','ref.system','ref.elapsed')

      runtimes[1:3] <- system.time(
        startPoints <- Jaatha.initialSearch(jaatha,nSim=200,nBlocksPerPar=4)
      )
      startPoints <- Jaatha.pickBestStartPoints(startPoints,best=2)

      runtimes[4:6] <- system.time(
        jaatha <- Jaatha.refineSearch(jaatha,startPoints,nSim=400,
                                    epsilon=.2,nFinalSim=200,
                                    halfBlockSize=0.05,weight=.9,
                                    nMaxStep=200)
      )
      estimates <- Jaatha.printLikelihoods(jaatha)[1,-(1:2)]
      sink(NULL)
      return(c(runtimes, estimates))
  }

  print(results)

  estimates <- results[, -(1:6)]
  runtimes <- results[, 1:6]

  write.table(estimates, file=paste(folder, "estimates.txt", sep=""), row.names=F)
  write.table(par.grid,  file=paste(folder, "true_values.txt", sep=""), row.names=F)
  write.table(runtimes,  file=paste(folder, "runtimes.txt", sep=""), row.names=F)
}

evalTest <- function(test.results){
  true.values <- test.results$RealPars
  estimates <- test.results$Estimates
  return(calcRSSE(estimates,true.values))
}

plotTrueVsEstimated <- function(true.values,estimates,parameter){
  plot(true.values[,parameter],estimates[,parameter])
  abline(0,1)
}

calcRSSE <- function(estimates,true.values){
  rsse <- sqrt(sum(estimates - true.values)^2)
  return(rsse)
}

calcAvgRunTime <- function(test.results){
  run.times <- apply(test.results$RunTimes[,-c(3,8)],1,sum)
  return(mean(run.times))
}

CreateParGrid <- function(dm,n.points){
  par.ranges <- jaatha:::dm.getParRanges(dm)
  #print(par.ranges)
  n.dim <- length(jaatha:::dm.getParameters(dm))
  n.runs <- n.points^n.dim
  par.grid <- matrix(0,n.runs,n.dim)
  
  for (i in 1:n.dim){
    pars.current.dim <- seq(par.ranges[i,1],par.ranges[i,2],length=n.points+2)[-c(1,n.points+2)]
    par.grid[,i] <- rep(pars.current.dim,each=n.points^(i-1))
  }

  return(par.grid)
}


dm.getParRanges <- function(dm,inklExtTheta=T){
    parMask <- !is.na(dm@features$parameter)
    parRanges <- cbind(lower=dm@features$lowerRange[parMask],upper=dm@features$upperRange[parMask])
    parRanges <- parRanges[sort.list(dm@features$parameter[parMask]),]
    if (!inklExtTheta) parRanges <- parRanges[1:dm.getNPar(dm),]
    return( parRanges )
}



#if (file.exists("testResults.save")) { 
#  load("testResults.save")
#} else { 
#  test.results <- list() 
#}

#if (is.null(test.results[[version]])) test.results[[version]] <- list()

#Test a simple theta/tau model
dm <- dm.createDemographicModel(c(24,25), 100)
dm <- dm.addSpeciationEvent(dm,0.001,5)
dm <- dm.addMutation(dm,1,20)
dm <- dm.addRecombination(dm,fixed=20)
dm <- dm.addSymmetricMigration(dm,fixed=.5)
runTest(dm, 5, model="tt")

#Test a model with 4 parameters
dm.mg <- dm.createDemographicModel(c(20,25), 100)
dm.mg <- dm.addSpeciationEvent(dm.mg,.001,5)
dm.mg <- dm.addMutation(dm.mg,1,20)
dm.mg <- dm.addSymmetricMigration(dm.mg, .1, 5)
dm.mg <- dm.addGrowth(dm.mg, .1, 5, population=2)
dm.mg <- dm.addRecombination(dm.mg, fixed=20)
runTest(dm.mg, 3, model="mg")
