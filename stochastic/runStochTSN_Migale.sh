#! /bin/bash

#
# This script submits jobs for the stochastic simulations.
#

# Compile the simulation script
gcc stochTSN_migale.c -lm -o stochTSN

# Path to the current folder
export THEPATH=DriveSpace/

# Number of replicates for each set of parameters
export NREPS=1

# Create execution files
# Make the execution files executable
# and submit them
for s in 0.4500000 0.4602273 0.4704545 0.4806818 0.4909091 0.5011364 0.5113636 0.5215909 0.5318182 0.5420455 0.5522727 0.5625000 0.5727273 0.5829545 0.5931818 0.6034091 0.6136364 0.6238636 0.6340909 0.6443182 0.6545455 0.6647727 0.6750000 0.6852273  0.6954545 0.7056818 0.7159091 0.7261364 0.7363636 0.7465909 0.7568182 0.7670455 0.7772727 0.7875000 0.7977273 0.8079545 0.8181818 0.8284091 0.8386364 0.8488636 0.8590909 0.8693182 0.8795455 0.8897727 0.9000000
do
  for mig in 0.1 0.5 1
  do
      # Create execution file
      echo -e "#!/bin/bash\n./${THEPATH}stochTSN 1000 ${s} ${mig} 10000 ${NREPS} > ${THEPATH}data/TSN_s-${s}_mig-${mig}.csv" > simTSN_s-${s}_mig-${mig}.sh
      # Make the file executable
      chmod +x simTSN_s-${s}_mig-${mig}.sh
      # Submit the file
      qsub -q short.q simTSN_s-${s}_mig-${mig}.sh
  done
done
