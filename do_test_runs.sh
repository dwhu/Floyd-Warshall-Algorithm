#!/bin/sh
##
nt=1
while [ $nt -le 10 ]
do
   MPI_WORKERS=$(($nt*$nt))
   echo "Running with "$MPI_WORKERS" workers"
   mpiexec -np $MPI_WORKERS floydWarshall A_3360.dat
   nt=`expr $nt + 1`
done