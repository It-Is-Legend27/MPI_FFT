#!/bin/bash

# Angel Badillo, Chad Callender
# Serial 1024 Frontera script

#SBATCH -J myserialjob_1024           # Job name
#SBATCH -o myserialjob_1024.o%j       # Name of stdout output file
#SBATCH -e myserialjob_1024.e%j       # Name of stderr error file
#SBATCH -p small           # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 00:2:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A ASC23018       # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=username@tacc.utexas.edu

# Any other commands must follow all #SBATCH directives...
module list
pwd
date

# Compile serial code
gcc Badillo_Callender_SeqVer_1024.c -lm -o serial_1024

# Launch serial code...
# loop runs code three times
for value in {1..3}
do
    ./serial_1024        # Do not use ibrun or any other MPI launcher
    echo ''              # newline
done