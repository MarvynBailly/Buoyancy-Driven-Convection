#!/bin/bash -l
#PBS -A UMCP0020
#PBS -N run_test
#PBS -o Logs/run_test.log
#PBS -j oe
#PBS -l walltime=05:00
#PBS -q casper
#PBS -l select=1:ncpus=1:ngpus=1:mem=8GB
#PBS -l gpu_type=v100
#PBS -M mbailly@umd.edu
#PBS -m abe

### Clear and load all the modules needed
module --force purge
module load ncarenv
module load cuda
module load julia/1.10.2

export TMPDIR=/glade/derecho/scratch/$USER/temp
proj_dir=/glade/u/home/mbailly/cuda_test3 #buoyancy_and_wind-driven_convection
mkdir -p $TMPDIR

### Run simulation
proj_dir=/glade/u/home/$USER/buoyancy_and_wind-driven_convection
time julia --project=$proj_dir simulations/buoyancy_and_wind-driven_convection.jl
###2>&1 | tee Logs/{0}.out

### Overwrite previous log file
LOG=$proj_dir/Logs/test.log
if [ -f "$LOG" ]; then
    rm -f $LOG
fi
mv $proj_dir/Logs/run_test.log $LOG

qstat -f $PBS_JOBID >> Logs/test.log

