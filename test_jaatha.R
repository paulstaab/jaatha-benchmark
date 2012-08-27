#!/usr/bin/Rscript --vanilla

# Tests the current installed version of Jaatha using
# a simple theta/tau model

library(jaatha)

version <- as.character(packageVersion("jaatha"))

cat("Testing version",version,"of Jaatha\n")

runTest <- function(dm, n.points=5, seed=12542){
  set.seed(seed)

  par.grid <- CreateParGrid(dm, n.points)

  cat("Simulating data...\n")
  sumStats <- dm.simSumStats(dm,par.grid)

  n <- dim(par.grid)[1]

  testResults <- list()
  testResults[['RealPars']] <- par.grid

  estimates <- matrix(0,n,dim(par.grid)[2])
  runtimes <- matrix(0,n,10)
  colnames(runtimes) <-
    c('init.user.self','init.sys.self','init.elapsed','init.user.child','init.sys.child',
      'ref.user.self','ref.sys.self','ref.elapsed','ref.user.child','ref.sys.child')

  for (i in 1:n){
      cat("----------------------------------------------------------------------\n")
      cat("Run",i,"of",n,"\n")
      cat("Real parameters:",par.grid[i,],"\n")
      cat("----------------------------------------------------------------------\n")
      jaatha <- Jaatha.initialize(dm, sumStats[i,])

      runtimes[i,1:5] <- system.time(
        startPoints <- Jaatha.initialSearch(jaatha,nSim=400,nBlocksPerPar=4)
      )
      startPoints <- Jaatha.pickBestStartPoints(startPoints,best=2)

      runtimes[i,6:10] <- system.time(
        jaatha <- Jaatha.refineSearch(jaatha,startPoints,nSim=400,
                                    epsilon=.2,nFinalSim=200,
                                    halfBlockSize=0.05,weight=.9,
                                    nMaxStep=200)
      )
      estimates[i,] <- Jaatha.printLikelihoods(jaatha)[1,-(1:2)]
  }

  testResults[['Estimates']] <- estimates
  testResults[['RunTimes']] <- runtimes

  return(testResults)
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





if (file.exists("testResults.save")) { 
  load("testResults.save")
} else { 
  test.results <- list() 
}

if (is.null(test.results[[version]])) test.results[[version]] <- list()

#Test a simple theta/tau model
dm <- dm.createDemographicModel(sampleSizes=c(24,25),nLoci=100,seqLength=1000)
dm <- dm.addSpeciationEvent(dm,0.001,5)
dm <- dm.addMutation(dm,1,20)
dm <- dm.addRecombination(dm,fixed=20)
dm <- dm.addSymmetricMigration(dm,fixed=.5)
test.results[[version]][['TauTheta']]  <- runTest(dm, 8)

#Test a model with 4 parameters
dm.mg <- dm.createDemographicModel(c(20,25), 100)
dm.mg <- dm.addSpeciationEvent(dm.mg,.001,5)
dm.mg <- dm.addMutation(dm.mg,1,20)
dm.mg <- dm.addSymmetricMigration(dm.mg, .1, 5)
dm.mg <- dm.addGrowth(dm.mg, .1, 5, population=2)
dm.mg <- dm.addRecombination(dm.mg, fixed=20)
test.results[[version]][['MigGrowth']]  <- runTest(dm, 4)

save(test.results, file="testResults.save")
