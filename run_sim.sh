#!/bin/bash -l

### Job Name
#PBS -N run_2D2-20-nosl
### Project Code Allocation
#PBS -A UMCP0020
### Resources
#PBS -l select=1:ncpus=16:mem=32GB:ngpus=1
### type of GPU
#PBS -l gpu_type=v100
### Run Time
#PBS -l walltime=03:00:00
#PBS -l job_priority=economy
### To the casper queue
###PBS -q gpudev
#PBS -q casper
###PBS -q develop
### Log file
### output
#PBS -o Logs/run_2D2-20-nosl.out
### error
#PBS -e Logs/run_2D2-20-nosl.err
### email 
#PBS -M mbailly@umd.edu
#PBS -m abe

### Module load
module --force purge
module --ignore-cache load ncarenv/23.10 gcc ncarcompilers netcdf
module --ignore-cache load cuda
module --ignore-cache load julia/1.9

### Run simulation
proj_dir=/glade/u/home/$USER/Buoyancy-Driven-Convection


julia --pkgimages=no --project=. Simulations/buoyancy_and_wind-driven_convection.jl --SL="no" --nTf=20 sim2D2-20-nosl

### Overwrite previous log file
#LOG=$proj_dir/Logs/run2D2.log
#if [ -f "$LOG" ]; then
#    rm -f $LOG
#fi
#mv $proj_dir/Logs/run_2D2.log $LOG
#
#qstat -f $PBS_JOBID >> Logs/run_2D2.log
