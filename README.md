A Benchmark Suite For Jaatha
===========

This is a benchmark for different version of 
[jaatha](https://github.com/paulstaab/jaatha). 
It test the performance  if Jaatha for a number of demographic models.

## Executing the benchmarks
The benchmarks are aimed to be executed on a Linux workstation using
```bash
./R/test_jaatha.R <path/to/jaatha_package.tar.gz>
```
If you want to execute the test on your own workstation, you need to adapt
the settings at the beginning of `test_jaatha.R`. It required to install
the 
[jaatha-benchmark-package](https://github.com/paulstaab/jaatha-benchmark-package)
first.

## Evaluting the results
The raw results are in the `runs` folder. For a quick evalution, call
```bash
./R/evaluate_results.R
```