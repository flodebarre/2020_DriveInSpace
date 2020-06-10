// To compile:
// gcc stochTSN.c -lm -o test

// To run:
// ./test 1000 0.575 0.1 10000 1 > simtext.csv
//         K     s   mig tfinal nrep

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h> // to get system time

// Simulation parameters
#define DTIME 0.05 // Time interval at which population is updated (has to be <1)

#define MAXSTEPS 10000000000 // Maximum number of simulation steps (in case a loop goes wrong)
#define STEPPRINTINTERVAL 10 // Print result every ... steps

#define NSITES 100 // Number of sites, i.e. discrete spatial locations

// Define variables
int nD[NSITES]; // number of DD individuals at each site
// NB: pop size is constant, so number of OO individuals is K-nD
int newnD[NSITES]; // temporary number of DD
double fecD; // Production of DD individuals
double fecO; // Production of OO individuals
int nLocal; // Temporary number
int nLeft; // Temporary number (migration)
int nRight; // Temporary number (migration)
int nStay; // Temporary number (migration)

int totDrive; // Total number of Drive alleles in the population

int timestep; // Time Step
double currenttime; // Time value
double deltaTime; // Time interval, if death occurs within DTIME

int isite; // Counter
int iindiv; // Counter
double randnum; // Random number

// Functions used in the script
void printPopulationState(void); // Print current state of the population
void oneStep(void); // Do one step of the simulation and update the population accordingly
void migrate(void); // substep: Migration

int g_argc; // Global variable for command line parameter entry
char **g_argv; // Global variable for command line parameter entry

// Introduce parameter names (values entered when running the script)
double K; // Carrying capacity
double s; // Zygote survival cost induced by the drive (0<= s <= 1)
double mig; // Emigration probability
double MAXTIME; // Maximum time until which simulation is run
int NREPS; // Number of simulation replicates

/* Function to process the parameters, give them names */
void process_command_line(void)
{
  // Extract parameters from the command line
  /// Model parameters
  K = (double) strtod(g_argv[1], (char **) NULL); // K, carrying capacity
  s = (double) strtod(g_argv[2], (char **) NULL); // Zygote survival cost induced by the drive (0<= s <=1)
  mig = (double) strtod(g_argv[3], (char **) NULL); // Emigration probability
  /// Simulation parameters
  MAXTIME = (double) strtod(g_argv[4], (char **) NULL); // Maximal time value of the simulation
  NREPS = (int) strtod(g_argv[5], (char **) NULL); // Number of replicates
}

int main(int argc, char **argv)
{
  // Process command line arguments
  g_argc = argc;
  g_argv = argv;
  process_command_line();
  printf("K = %.0f, s = %.4f, mig = %.4f, MAXTIME = %.2f, NREPS = %d\n", K, s, mig, MAXTIME, NREPS);

  // Set random seed to ensure reproducibility
  unsigned long seed=98756412; // to be changed
  srand48(seed);

  int irep = 0; // REPLICATE number

  for(irep = 0; irep<((int)NREPS); irep++){

    // Initialize time
    currenttime = 0.0;

    // Initialize the population
      // WT on the left-hand side, at carrying capacity K
    for(isite=0; isite< (int)floor(NSITES/2); isite++){
      nD[isite] = (int) 0;
    }
      // Drive on the right-hand side, at K too
      // NB: this is not the carrying capacity of an all-drive population when s>0
    for(isite = (int)floor(NSITES/2); isite < NSITES; isite++){
      nD[isite] = (int) floor(K);
    }

//    // If only one replicate, print population state
//    if(NREPS<2){ // If one sim, print output at regular time intervals
//      printPopulationState();
//    }

    // Compute total population size and total drive
    totDrive = 0;
    for(isite=0; isite<NSITES; isite++){
      totDrive += nD[isite];
    }

    for(timestep = 0; timestep<((int)MAXSTEPS); timestep++){
      // Check that the population is not extinct, that the drive is still present
      // and that time is not over
      if((NSITES*K - totDrive)<1 || totDrive<1 || currenttime >= MAXTIME){
          break;
      }
      // Do one step of the simulation
      oneStep();

//      // If only one replicate, print population state
//      if(NREPS<2){ // If one sim, print output at regular time intervals
//        if(1.0*timestep - STEPPRINTINTERVAL * floor(timestep/STEPPRINTINTERVAL) < 1.0){
//          // if round number of timestep in base STEPPRINTINTERVAL, print ouput
//          printPopulationState();
//        }
//      }
    }

    // Simulation finished -> save final outcome
    printPopulationState();
  } //end rep
  
  return(0);

} // end main Functions



/*----------------------------------------------------------------------------*/
/* ONE SIMULATION STEP - ZYGOTE CONVERSION */
void oneStep(void){
  // Update time
  currenttime += (1.0)*DTIME;
  
  //-------------------------------------------
  // 1) Death and birth
  // Update each site one at a time
  for(isite=0; isite<NSITES; isite++){
    // Reset newpop
    newnD[isite] = 0;

    // Compute fecs
    fecD = (double) 1.0 * (1.0 - s) * (nD[isite] * (2.0 * (K - nD[isite]) + nD[isite]))/(1.0 * K * K);
    fecO = (double) 1.0 * (K - nD[isite]) * (K - nD[isite])/(1.0 * K * K);

    // For each individual, draw its type
    for(iindiv=0; iindiv<(int) floor(K); iindiv++){
      randnum = drand48();
      if(randnum < (double)(fecD / (fecD+fecO))){
        // New individual is Drive
        newnD[isite] +=1;
      }
      // Otherwise, new individual is WT (no need to track them since total size is fixed =K)
    }
  }
  
  //--------------------------------------------
  // 2) Migration

  // Reset pop
  for(isite=0; isite<NSITES; isite++){
    nD[isite] = 0;
  }

  // Each site, excluding edges
  for(isite=1; isite<(NSITES-1); isite++){
    // Drive individuals
    nLocal = newnD[isite];
    migrate();
    nD[isite-1] += nLeft;
    nD[isite+1] += nRight;
    nD[isite] += nStay;
  }
  // Edges
  // Site 0 -> to the left in fact stay home
    // Drive individuals
    nLocal = newnD[0];
    migrate();
    nD[0] += nLeft + nStay;
    nD[1] += nRight;

  // Site (NSITES-1) -> to the right in fact stay home
    // Drive individuals
    nLocal = newnD[NSITES-1];
    migrate();
    nD[NSITES-1] += nRight + nStay;
    nD[NSITES-2] += nLeft;

  //--------------------------------------------
  // Updates
  totDrive = 0;

  for(isite=0; isite<NSITES; isite++){
    totDrive += nD[isite];
  }
}

/*----------------------------------------------------------------------------*/
/* SUBSTEP: MIGRATION */
void migrate(void){
  // For a specified site and a specified type, nLocal
  
  // Reset output variables
  nLeft = 0;
  nRight = 0;
  nStay = 0;
  
  // Here we pick a fixed proportion of individuals who will migrate
  // Just pick whether floor or ceiling
  randnum = drand48();
  if(randnum < 0.5){
    nStay = (int) floor((1.0-mig) * nLocal);
  }else{
    nStay = (int) ceil((1.0-mig) * nLocal);
  }
  
  randnum = drand48();
  if(randnum < 0.5){
    nLeft = (int) floor((1.0 * nLocal - nStay)/2.0);
  }else{
    nLeft = (int) ceil((1.0 * nLocal - nStay)/2.0);
  }

  nRight = nLocal - nStay - nLeft;
}

/*----------------------------------------------------------------------------*/
/* IF ONE REP: PRINT POPULATION STATE */
/* PRINT STATE OF THE POPULATION */
void printPopulationState(void){
  // Print current time
  printf("%09.3f", currenttime);
  // Print number of individuals of each genotype
  for (isite=0; isite<NSITES; isite++) {
    printf(", %d", nD[isite]);
  }
  printf("\n");
}

