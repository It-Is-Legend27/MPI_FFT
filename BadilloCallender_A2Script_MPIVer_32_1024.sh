#!/bin/bash

# Angel Badillo, Chad Callender
# Parallel (32, 1024) Frontera script

#SBATCH -J myparalleljob_32_1024           # Job name
#SBATCH -o myparalleljob_32_1024.o%j       # Name of stdout output file
#SBATCH -e myparalleljob_32_1024.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 4               # Total # of nodes 
#SBATCH -n 32              # Total # of mpi tasks
#SBATCH -t 00:2:00         # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A ASC23018       # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=username@tacc.utexas.edu

# Any other commands must follow all #SBATCH directives...
module list
pwd
date

# Compile parallel code
mpicc Badillo_Callender_MPIVer_1024.c -lm -o parallel_32_1024

# Launch parallel code...
# loop runs code three times
for value in {1..3}
do
    ibrun ./parallel_32_1024       # Use ibrun instead of mpirun or mpiexec
    echo ''              # newline
done
