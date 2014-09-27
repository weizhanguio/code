#!/bin/bash
#PBS -A nn2849k
#PBS -N job
#PBS -e err

#PBS -l walltime=1:00:00
#PBS -l select=16:ncpus=16:mpiprocs=1:ompthreads=16
#### select=8:ncpus=16:mpiprocs=2:ompthreads=8
 
# Tips to improving performance
#
# 1. Adjusting MPI_BUFS_PER_PROC and MPI_BUFS_PER_HOST.
#
# Use the "-stats" option to mpiexec_mpt to get additional information in the
# output file. Included in that information is the number of retries for
# allocating MPI buffers. After you have executed your program with the "-stats"
# option, you can see this by typing something similar to:
#
# $ cat my_mpi_job.o55809 | grep retries | grep -v " 0 retries"
#
# You can then increase the values of MPI_BUFS_PER_PROC (default 32) and
# MPI_BUFS_PER_HOST (default 96) until the number of retries is sufficiently
# low, e.g. by uncommenting these lines:
#
# export MPI_BUFS_PER_PROC=256
# export MPI_BUFS_PER_HOST=1024
#
# See "man mpi" for more information.
#
#
# 2. Adjusting MPI_BUFFER_MAX
#
# For some codes it gives a significant increase in performance to specify a
# value for MPI_BUFFER_MAX. According to "man mpi" this value "Specifies a
# minimum message size, in bytes, for which the message will be considered a
# candidate for single-copy transfer." The value of MPI_BUFFER_MAX varies from
# program to program, but typical values are between 2048 and 32768. You can
# therefore test if this improves the performance of your program by executing
# it like this:
#
# export MPI_BUFFER_MAX=2048
# time -p mpiexec_mpt -n 1024 ./myprog
#
# See "man mpi" for more information.
 
module load mpt
#module load intelcomp/14.0.1 
module load intelcomp/12.0.5.220
#module load scalasca/1.4.2
cd $PBS_O_WORKDIR
#export MPI_DSM_CPULIST=0,15:allhosts
#export MPI_DSM_CPULIST=0,1,12,13:allhosts
#export MPI_DSM_CPULIST=0,1,2,3,10,11,12,13:allhosts
#export MPI_DSM_CPULIST=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15:allhosts

#export MPI_DSM_VERBOSE=1 
export MPI_BUFFER_MAX=32768
#mpiexec_mpt -n 4  omplace -nt 16    ./app  100  100 100 400
mpiexec_mpt -n  16   omplace -nt 16  ./app   100 100 100  3000     2  2  4
#mpirun -np 2 omplace  -vv ./app 0.1 2000
#export KMP_AFFINITY disabled
#export OMP_NUM_THREADS=2
#omplace -vv -c 0,15 ./app   60  60 60 2000    1  1  1
#export OMP_NUM_THREADS=4
#omplace -vv -c 0,1,11,12  ./app  60  60 60 2000    1  1  1
#export OMP_NUM_THREADS=8
#omplace  -c 0,1,2,3,10,11,12,13  -vv ./app  60  60 60 2000    1  1  1
#export OMP_NUM_THREADS=16
#omplace -vv -c 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15  ./app   100  100 100 3000    1  1  1


