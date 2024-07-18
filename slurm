#!/bin/bash
#SBATCH --nodes=16           # number of nodes
#SBATCH --ntasks-per-node=128   # number of processors per node
#SBATCH -t 24:00:00      # run for 3 hr max
#SBATCH --mem=499gb
#SBATCH --tmp=100gb
#SBATCH -e stderr              # Name of output file for error messages
#SBATCH -o stdout              # Name of output file for standard output
#SBATCH --partition=crfel
#SBATCH --job-name=NDW_3D
##SBATCH --nodelist=acn[151,157-158,200,206-207,216,218,220-221,227-230,236]
##SBATCH --nodelist=acn[157-158,206-207,216,218,220-221]

module purge

module load cmake/3.26.3
module load ompi/4.1.5/gnu-11.3.0

cd /home/suo-yang/ramac106/run/NDW_2D/H2_Salinas

#make realclean TPLrealclean
#make TPL

make clean
make -j2048 2>&1 | tee surya.build.log

# ignore the slurm/pmix machinery
export OMPI_MCA_plm=rsh
export OMPI_MCA_sec=^munge

# make sure to pass LD_LIBRARY_PATH over ssh
mpirun -x LD_LIBRARY_PATH -np 2048 ./PeleC2d.gnu.MPI.ex inp.3d 2>&1 | tee surya.run.log
