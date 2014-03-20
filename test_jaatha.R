#!/usr/bin/Rscript --vanilla

arg <- commandArgs()
path <- '/scratch/paul/test_jaatha/lib'

if (!is.na(arg[2])) {
  stopifnot(file.exists(arg[2]))
  jaatha.package <- arg[2]
} else {
  jaatha.package <- 'packages/jaatha_2.97.2.tar.gz'
}

dir.create(path, recursive=TRUE)
install.packages(jaatha.package, dependencies=TRUE, lib=path) 
library(jaatha, lib.loc=path)
library(testJaatha)

version <- as.character(packageVersion("jaatha"))

# Test a simple theta/tau/migration model
dm <- dm.createDemographicModel(c(20,25), 75)
dm <- dm.addSpeciationEvent(dm, 0.01, 5)
dm <- dm.addMutation(dm, 1, 10)
dm <- dm.addRecombination(dm, fixed=5)
dm <- dm.addSymmetricMigration(dm, .1, 2)
testJaatha:::testJaatha(dm, 3, 2, seed=12579, smoothing=FALSE, cores=c(8,4),
                        folder=paste('runs', version, 'tt.classic', sep='/'))
testJaatha:::testJaatha(dm, 3, 2, seed=12579, smoothing=TRUE, cores=c(8,4),
                        folder=paste('runs', version, 'tt.smooth', sep='/'))


# Test the fpc sum.stat
dm.fpc <- dm.createDemographicModel(c(20,25), 100)
dm.fpc <- dm.addSpeciationEvent(dm.fpc, 0.01, 5)
dm.fpc <- dm.addMutation(dm.fpc, 1, 10)
dm.fpc <- dm.addRecombination(dm.fpc, 1, 10)
dm.fpc <- dm.addSymmetricMigration(dm.fpc, .1, 2)
testJaatha:::testJaatha(dm.fpc, 2, 2, seed=124578, smoothing=TRUE, cores=c(8, 4),
                        folder=paste('runs', version, 'fpc', sep='/'),
                        fpc=TRUE)


# Test a model with 5 parameters
dm.mg <- dm.createDemographicModel(c(20,25), 100)
dm.mg <- dm.addSpeciationEvent(dm.mg,.001,5)
dm.mg <- dm.addMutation(dm.mg,1,20)
dm.mg <- dm.addMigration(dm.mg, .1, 5, pop.from=1, pop.to=2,
                         new.par.name="m12")
dm.mg <- dm.addMigration(dm.mg, .1, 5, pop.from=2, pop.to=1, 
                         new.par.name="m21")
dm.mg <- dm.addGrowth(dm.mg, .1, 5, population=2)
dm.mg <- dm.addRecombination(dm.mg, fixed=20)
testJaatha:::testJaatha(dm.mg, 2, 2, seed=24680, smoothing=FALSE, cores=(32,1),
                        folder=paste('runs', version, 'mg.classic', sep='/'))
testJaatha:::testJaatha(dm.mg, 2, 2, seed=24680, smoothing=TRUE, cores=(32,1),
                        folder=paste('runs', version, 'mg.smooth', sep='/'))
