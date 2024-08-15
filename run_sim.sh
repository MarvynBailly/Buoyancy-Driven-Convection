#!/bin/bash -l

### Job Name
#PBS -N run_2D2
### Project Code Allocation
#PBS -A UMCP0020
### Resources
#PBS -l select=1:ncpus=16:mem=32GB:ngpus=1
### Run Time
#PBS -l walltime=03:00:00
### To the casper queue
###PBS -q gpudev
#PBS -q casper
### output
#PBS -o Logs/run_2D2.out
### error
#PBS -e Logs/run_2D2-20.err
### type of GPU
#PBS -l gpu_type=v100
### email 
#PBS -M mbailly@umd.edu
#PBS -m abe


module --force purge
module --ignore-cache load ncarenv/23.10 gcc ncarcompilers netcdf
module --ignore-cache load cuda
module --ignore-cache load julia/1.9
### file to run

julia --pkgimages=no --project=. Simulations/buoyancy_and_wind-driven_convection.jl --outdir="./Data/" 2D2
