#!/bin/bash
#SBATCH -J Serial_cpu_job
#SBATCH -p batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --comment inhouse #Refer to the Table of SBATCH option name per application
 
 
./a.out 
 
exit 0
