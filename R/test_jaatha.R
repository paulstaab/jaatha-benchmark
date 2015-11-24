#!/usr/bin/Rscript --no-save

library(jaatha)
library(coala)
library(testJaatha)
library(devtools) # to keep it in packrat

version <- packageVersion("jaatha")


# --- Test a tree population model ---------------------------------------------
model_3pop <- coal_model(c(15, 15, 15), 200) +
  feat_pop_merge(par_range("tau21", 0.01, 5), 2, 1) +
  feat_pop_merge(par_range("tau31", 0.01, 5), 3, 1) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_recombination(5) +
  feat_migration(par_range("m", .1, 2), symmetric = TRUE) +
  sumstat_jsfs(population = 1:3)

testJaatha(model_3pop, 2, 2, seed = 12542, cores = c(16, 2),
           folder = file.path('runs', version, '3pop'))


# --- Test a simple theta/tau/migration model ---------------------------------
model <- coal_model(c(20, 25), 75) +
  feat_pop_merge(par_range("tau", 0.01, 5), 2, 1) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_recombination(5) +
  feat_migration(par_range("m", .1, 2), symmetric = TRUE) +
  sumstat_jsfs()

set.seed(124)
test_data <- createTestData(model, 3, 2, cores=32)
testJaatha(model, test_data = test_data, seed = 125, cores=c(16, 2),
           folder=file.path('runs', version, 'tt_is'))
testJaatha(model, test_data = test_data, seed = 126, cores=c(16, 2),
           init_method = "zoom-in",
           folder=file.path('runs', version, 'tt_zoom'))
testJaatha(model, test_data = test_data, seed = 127, cores=c(16, 2),
           init_method = "middle",
           folder=file.path('runs', version, 'tt_middle'))


# --- Test the fpc sum.stat ---------------------------------------------------
model_fpc <- coal_model(c(20,25), 100) +
  feat_pop_merge(par_range("tau", 0.01, 5), 2, 1) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_recombination(par_range("rho", 1, 10)) +
  feat_migration(par_range("m", .1, 2), symmetric = TRUE) +
  sumstat_jsfs() +
  sumstat_four_gamete("fgc1", 1) +
  sumstat_four_gamete("fgc2", 2)

testJaatha(model_fpc, 2, 2, seed=2537, cores=c(16, 2),
           folder=file.path('runs', version, 'fpc'))


# --- Test a finite sites model -----------------------------------------------
model_fs <- coal_model(c(15, 20, 2), 50) +
  feat_pop_merge(par_range("tau", 0.01, 5), 2, 1) +
  feat_mutation(par_range("theta", 1, 10), model = "HKY", 
                base_frequencies = c(0.2, 0.2, 0.3, 0.3),
                tstv_ratio = 2) +
  feat_recombination(1) +
  feat_migration(par_range("m", .1, 2), symmetric = TRUE) +
  feat_outgroup(3) +
  feat_pop_merge(par_expr("2*tau"), 3, 1) +
  sumstat_jsfs()

testJaatha(model_fs, 2, 4, seed=124578, cores=c(16, 2),
           folder=file.path('runs', version, 'fs'))


# --- Test a 7 par sites model ------------------------------------------------
model_6par <- coal_model(c(15, 20), 500, 1000) +
  feat_mutation(par_range('theta', .1, 10)) +
  feat_migration(par_range('m12', 0.01, 5), 1, 2) +
  feat_migration(par_range('m21', 0.01, 5), 2, 1) +
  feat_size_change(par_range('q', 0.01, 10), population = 2, time = 0) +
  par_range("s2", 0.001, 2) +
  feat_growth(par_expr(log(1)/tau), population = 1, time = 0) +
  feat_growth(par_expr(log(q/s2)/tau), population = 2, time = 0) +
  feat_size_change(par_expr(1 + s2), population = 1, time = par_expr(tau)) +
  feat_pop_merge(par_range('tau', 0.001, 5), 2, 1) +
  feat_recombination(par_const(5)) +
  sumstat_jsfs()

testJaatha(model_6par, 2, 1, seed = 24680, cores = c(16, 2),
           folder = file.path('runs', version, '6par'))
#testJaatha(model_6par, 2, 1, seed=24680, smoothing=TRUE, cores=c(16, 2),
#           folder=file.path('runs', version, '6par_smoothing'))
