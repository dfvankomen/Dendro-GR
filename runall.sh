#!/bin/bash

trap 'exit 120' INT

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_BLOSC_double_Level3.par.toml | tee bssnq1_BLOSC_double_Level3_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_BLOSC_double_Level6.par.toml | tee bssnq1_BLOSC_double_Level6_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_BLOSC_double_Level9.par.toml | tee bssnq1_BLOSC_double_Level9_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_double_accuracy0.par.toml | tee bssnq1_ZFP_double_accuracy0_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_double_accuracy1.par.toml | tee bssnq1_ZFP_double_accuracy1_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_double_accuracy2.par.toml | tee bssnq1_ZFP_double_accuracy2_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_double_accuracy3.par.toml | tee bssnq1_ZFP_double_accuracy3_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_double_precision4.par.toml | tee bssnq1_ZFP_double_precision4_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_double_precision5.par.toml | tee bssnq1_ZFP_double_precision5_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_double_precision6.par.toml | tee bssnq1_ZFP_double_precision6_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_double_precision7.par.toml | tee bssnq1_ZFP_double_precision7_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_double_precision8.par.toml | tee bssnq1_ZFP_double_precision8_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_BLOSC_float_Level3.par.toml | tee bssnq1_BLOSC_float_Level3_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_BLOSC_float_Level6.par.toml | tee bssnq1_BLOSC_float_Level6_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_BLOSC_float_Level9.par.toml | tee bssnq1_BLOSC_float_Level9_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_float_accuracy0.par.toml | tee bssnq1_ZFP_float_accuracy0_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_float_accuracy1.par.toml | tee bssnq1_ZFP_float_accuracy1_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_float_accuracy2.par.toml | tee bssnq1_ZFP_float_accuracy2_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_float_accuracy3.par.toml | tee bssnq1_ZFP_float_accuracy3_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_float_precision4.par.toml | tee bssnq1_ZFP_float_precision4_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_float_precision5.par.toml | tee bssnq1_ZFP_float_precision5_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_float_precision6.par.toml | tee bssnq1_ZFP_float_precision6_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_float_precision7.par.toml | tee bssnq1_ZFP_float_precision7_out.txt 

mpirun -np 12 ./bssnSolver ../../parstest/bssnq1_ZFP_float_precision8.par.toml | tee bssnq1_ZFP_float_precision8_out.txt 

