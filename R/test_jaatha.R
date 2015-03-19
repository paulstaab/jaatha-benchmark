#!/usr/bin/Rscript --no-save --no-site-file --no-init-file

library(jaatha)
library(coalsimr)
library(testJaatha)

version <- paste0(packageVersion("jaatha"), "_",
                  packageVersion("coalsimr"))

# Test a simple theta/tau/migration model
model <- coal_model(c(20,25), 75) +
  feat_pop_merge(par_range("tau", 0.01, 5), 2, 1) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_recombination(5) +
  feat_migration(par_range("m", .1, 2), symmetric = TRUE) +
  sumstat_jsfs()

testJaatha(model, 3, 2, seed=12579, smoothing=FALSE, cores=c(16, 2),
           folder=file.path('runs', version, 'tt.old'))
testJaatha(model, 3, 2, seed=12579, smoothing=TRUE, cores=c(16, 2),
           folder=file.path('runs', version, 'tt.sm'))

# Test the fpc sum.stat
model_fpc <- coal_model(c(20,25), 100) +
  feat_pop_merge(par_range("tau", 0.01, 5), 2, 1) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_recombination(par_range("rho", 1, 10)) +
  feat_migration(par_range("m", .1, 2), symmetric = TRUE) +
  sumstat_jsfs() +
  sumstat_four_gamete("fgc1", 1) +
  sumstat_four_gamete("fgc2", 2)
testJaatha(model_fpc, 3, 2, seed=2537, smoothing=FALSE, cores=c(16, 2),
           folder=file.path('runs', version, 'tt.old'))


# Test a finite sites model
# dm.fs <- dm.setMutationModel(dm.fpc, "HKY", c(0.2, 0.2, 0.3, 0.3), 2)
# dm.fs <- dm.addOutgroup(dm.fs, '2*t_split_1', 2)
# testJaatha:::testJaatha(dm.fs, 2, 3, seed=124578, smoothing=FALSE, cores=c(16, 2),
#                         folder=paste('runs', version, 'fs', sep='/'))

# Test a model with 5 parameters
# dm.gro <- dm.createDemographicModel(c(24,25), 100, 1000)
# dm.gro <- dm.addMutation(dm.gro, 1, 20, parameter="theta")
# dm.gro <- dm.addSpeciationEvent(dm.gro, 0.001, 20, time.point="tau")
# dm.gro <- dm.addSymmetricMigration(dm.gro, 0.005, 5, parameter="m")
# dm.gro <- dm.addRecombination(dm.gro, parameter=5)
# dm.gro <- dm.addSizeChange(dm.gro, 0.05, 10, population=2, at.time="0", parameter="q")
# dm.gro <- dm.addParameter(dm.gro, 0.05, 10, par.name="s1")
# dm.gro <- dm.addParameter(dm.gro, 0.05, 10, par.name="s2")
# dm.gro <- dm.addSizeChange(dm.gro, parameter="s1+s2", population=1, at.time="tau", min.size.factor=NA)
# dm.gro <- dm.addGrowth(dm.gro, parameter="log(1/s1)/tau", population=1)
# dm.gro <- dm.addGrowth(dm.gro, parameter="log(1/s2)/tau", population=2)
# testJaatha:::testJaatha(dm.gro, 2, 1, seed=24680, smoothing=FALSE, cores=c(16, 2),
#                         folder=paste('runs', version, 'gro.old', sep='/'))
#testJaatha:::testJaatha(dm.gro, 2, 1, seed=24680, smoothing=TRUE, cores=c(16, 2),
#                        folder=paste('runs', version, 'gro.sm', sep='/'))

