#!/bin/bash
#PBS -A nn2849k
#PBS -N app
#PBS -e err
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=16:mpiprocs=16
 
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
 export MPI_BUFFER_MAX=32768
# time -p mpiexec_mpt ./myprog
#
# See "man mpi" for more information.
 
module load mpt
module load intelcomp/12.0.5.220 
 
cd $PBS_O_WORKDIR
 

#export MPI_DSM_CPULIST=0,15:allhosts
#mpiexec_mpt -n  2 ./app    1000 10000
#mpiexec_mpt -n 2 ./app    2000 20000
#mpiexec_mpt -n 2 ./app    4000 40000


#export MPI_DSM_CPULIST=0,1,12,13:allhosts
#mpiexec_mpt -n  4 ./app    1000 10000
#mpiexec_mpt -n 4 ./app    2000 20000
#mpiexec_mpt -n 4 ./app    4000 40000


#export MPI_DSM_CPULIST=0,1,2,3,10,11,12,13:allhosts
#mpiexec_mpt -n  8 ./app    1000 10000
#mpiexec_mpt -n 8 ./app    2000 20000
#mpiexec_mpt -n 8 ./app    4000 40000

export MPI_DSM_CPULIST=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15:allhosts
#mpiexec_mpt -n  16 ./app    1000 10000
mpiexec_mpt -n 16 ./app    2000 20000
#mpiexec_mpt -n 16 ./app    4000 40000




