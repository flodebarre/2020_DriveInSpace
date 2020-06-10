#! /bin/bash

#
# This script submits jobs for the stochastic simulations.
#

# Compile the simulation script
gcc stochwave_migale.c -lm -o stochwave

# Path to the current folder
export THEPATH=DriveSpace/

# Number of replicates for each set of parameters
export NREPS=1

# Create execution files
# Make the execution files executable
# and submit them
for s in 0.4500000 0.4604167 0.4708333 0.4812500 0.4916667 0.5020833 0.5125000 0.5229167 0.5333333 0.5437500 0.5541667 0.5645833 0.5750000 0.5854167 0.5958333 0.6062500 0.6166667 0.6270833 0.6375000 0.6479167 0.6583333 0.6687500 0.6791667 0.6895833 0.7000000
do
  for r in 0.12 0.24 0.36 0.48 0.60 0.72 0.84 0.96 1.08 1.20 1.32 1.44 1.56 1.68 1.80 1.92 2.04 2.16 2.28 2.40 2.52 2.64 2.76 2.88 3.00
  do
    for mig in 0.1 1
    do
      # Create execution file
      echo -e "#!/bin/bash\n./${THEPATH}stochwave 1000 ${s} ${r} ${mig} 10000 ${NREPS} > ${THEPATH}data/stoch_s-${s}_r-${r}_mig-${mig}.csv" > sim_s-${s}_r-${r}_mig-${mig}.sh
      # Make the file executable
      chmod +x sim_s-${s}_r-${r}_mig-${mig}.sh
      # Submit the file
      qsub -q short.q sim_s-${s}_r-${r}_mig-${mig}.sh
    done
  done
done
