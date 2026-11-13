#!/bin/sh
#SBATCH --output=AlNbTiVZr-shearYZp04.sout
#SBATCH --job-name=MS_QE_Elastic_AlNbTiVZr-shearYZp04
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=All
#source /opt/intel/oneapi/setvars.sh
rm -rf pwscf*
node=1
threads=1
processors=$(nproc)
np=$(($node*$processors/$threads))

export OMP_NUM_THREADS=$threads
export MKL_DEBUG_CPU_TYPE=5
export MKL_CBWR=AUTO
export LD_LIBRARY_PATH=/opt/mpich-4.0.3/lib:/opt/intel/oneapi/mkl/latest/lib:$LD_LIBRARY_PATH
export PATH=/opt/mpich-4.0.3/bin:$PATH
mpiexec -np $np /opt/thermoPW-7-2_intel/bin/pw.x -in AlNbTiVZr-shearYZp04.in
rm -rf pwscf*
perl /opt/qe_perl/QEout_analysis.pl
perl /opt/qe_perl/QEout2data.pl
