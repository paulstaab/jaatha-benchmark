#!/usr/bin/Rscript --vanilla

# Tests the current installed version of Jaatha using
# a simple theta/tau model

library(jaatha)
library(foreach)
library(doMC)

result.folder <- "results"
version <- as.character(packageVersion("jaatha"))

cat("Testing version",version,"of Jaatha\n")

runTest <- function(dm, n.points=2, seed=12523, model){
  cat("Testing", model, "\n")
  set.seed(seed)

  par.grid <- createParGrid(dm, n.points)

  cat("Simulating data...\n")

  n <- dim(par.grid)[1]
  seeds <- sample.int(2^20, n)

  folder <- paste("results", version, model, sep="/")
  log.folder <- paste("logs", version, model, sep="/")
  dir.create(folder, recursive=T, showWarnings=F)
  dir.create(log.folder, recursive=T, showWarnings=F)

  registerDoMC(8)

  results <- foreach(i=1:n, .combine=rbind) %dopar% { 
      cat("Run",i,"of",n,"\n")
      sink(file=paste(log.folder, "/run_", i, ".txt", sep=""))
      set.seed(seeds[i])
      cat("----------------------------------------------------------------------\n")
      cat("Run",i,"of",n,"\n")
      cat("Real parameters:",par.grid[i,],"\n")
      cat("----------------------------------------------------------------------\n")
      jsfs <- dm.simSumStats(dm, par.grid[i, ])
      jaatha <- Jaatha.initialize(dm, jsfs=jsfs, cores=4)

      runtimes <- rep(0, 6)
      names(runtimes) <-
        c('init.user','init.system','init.elapsed','ref.user','ref.system','ref.elapsed')

      runtimes[1:3] <- system.time(
        jaatha <- Jaatha.initialSearch(jaatha, 100, 2)
      )

      runtimes[4:6] <- system.time(
        jaatha <- Jaatha.refinedSearch(jaatha, 2 , 100, 100, epsilon=1)
      )
      estimates <-  Jaatha.getLikelihoods(jaatha)[1,-(1:2)]
      sink(NULL)
      return(c(runtimes, estimates))
  }

  estimates <- results[, -(1:6)]
  runtimes <- results[, 1:6]

  write.table(estimates, file=paste(folder, "estimates.txt", sep="/"), row.names=F)
  write.table(par.grid,  file=paste(folder, "true_values.txt", sep="/"), row.names=F)
  write.table(runtimes,  file=paste(folder, "runtimes.txt", sep="/"), row.names=F)
}

createParGrid <- function(dm,n.points){
  par.ranges <- jaatha:::dm.getParRanges(dm)
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


#Test a simple theta/tau/migration model
dm <- dm.createDemographicModel(c(20,25), 100)
dm <- dm.addSpeciationEvent(dm, 0.01, 5)
dm <- dm.addMutation(dm, 1, 10)
dm <- dm.addRecombination(dm,fixed=20)
dm <- dm.addSymmetricMigration(dm, .1, 2)
runTest(dm, 3, model="tt", seed=13579)

#Test a model with 5 parameters
dm.mg <- dm.createDemographicModel(c(20,25), 100)
dm.mg <- dm.addSpeciationEvent(dm.mg,.001,5)
dm.mg <- dm.addMutation(dm.mg,1,20)
dm.mg <- dm.addMigration(dm.mg, .1, 5, pop.from=1, pop.to=2,
                         new.par.name="m12")
dm.mg <- dm.addMigration(dm.mg, .1, 5, pop.from=2, pop.to=1, 
                         new.par.name="m21")
dm.mg <- dm.addGrowth(dm.mg, .1, 5, population=2)
dm.mg <- dm.addRecombination(dm.mg, fixed=20)
#runTest(dm.mg, 2, model="mg", seed=24680)
