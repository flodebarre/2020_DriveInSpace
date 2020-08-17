// To compile:
// gcc stochWolbachia.c -lm -o test

// To run:
//////////////////// TO BE UPDATED
// ./test 0.3 0.1 1000 1.25 0.1 10000 1 > simtext.csv
//        fw  wh  K    r   mig tfinal nrep

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
int nW[NSITES]; // number of Wolbachia-carrying individuals at each site
int nS[NSITES]; // number of non-infected individuals at each site
int newnW[NSITES]; // temporary number of W
int newnS[NSITES]; // temporary number of S
int newFocal; // Whether the individual that we are considering survives
int newOffspring; // Number of offspring produced by the focal
double fecundity; // Expected number of offspring in the time interval
int nLocal; // Temporary number
int nLeft; // Temporary number (migration)
int nRight; // Temporary number (migration)
int nStay; // Temporary number (migration)

int totPopSize; // Total population size
int totW; // Total number of individuals with Wolbachia in the population
int totS; // Total number of uninfected individuals

int timestep; // Time Step
double currenttime; // Time value
double deltaTime; // Time interval, if death occurs within DTIME

int isite; // Counter
int iindiv; // Counter
double randnum; // Random number
double randDeathTime; // Random time

double pPois; // Temporary variable for Poisson sampling
double LPois; // Temporary variable for Poisson sampling

// Functions used in the script
void printPopulationState(void); // Print current state of the population
void oneStep(void); // Do one step of the simulation and update the population accordingly
void deathAndBirth(void); // substep: Death and Birth
void migrate(void); // substep: Migration

int g_argc; // Global variable for command line parameter entry
char **g_argv; // Global variable for command line parameter entry

// Introduce parameter names (values entered when running the script)
double fw; // Factor affecting fecundity of Wolbachia-carrying mothers
double wh; // Factor affecting survival of offspring of W mother and S father
double K; // Carrying capacity
double r; // Growth rate
double mig; // Emigration probability
double MAXTIME; // Maximum time until which simulation is run
int NREPS; // Number of simulation replicates

// Other parameters
double deathRate = 1.0; // Scaled parameter

/* Function to process the parameters, give them names */
void process_command_line(void)
{
  // Extract parameters from the command line
  /// Model parameters
  fw = (double) strtod(g_argv[1], (char **) NULL); // factor affecting fertility of Wolbachia-carrying mothers
  wh = (double) strtod(g_argv[2], (char **) NULL); // factor affecting survival of offspring of infected mothers and uninfected fathers
  K = (double) strtod(g_argv[3], (char **) NULL); // K, carrying capacity
  r = (double) strtod(g_argv[4], (char **) NULL); // Growth rate
  mig = (double) strtod(g_argv[5], (char **) NULL); // Emigration probability
  /// Simulation parameters
  MAXTIME = (double) strtod(g_argv[6], (char **) NULL); // Maximal time value of the simulation
  NREPS = (int) strtod(g_argv[7], (char **) NULL); // Number of replicates
}

int main(int argc, char **argv)
{
  // Process command line arguments
  g_argc = argc;
  g_argv = argv;
  process_command_line();
  printf("fw = %.4f, wh = %.4f, K = %.0f, r = %.2f, mig = %.4f, MAXTIME = %.2f, NREPS = %d\n", fw, wh, K, r, mig, MAXTIME, NREPS);

  // Random seed
  if(NREPS>1){ // If multiple replicates
    // Set random seed to ensure reproducibility
    unsigned long seed=98756412; // to be changed
    srand48(seed);
  }else{
    // or change it everytime
    srand48(time(NULL));
  }

  int irep = 0; // REPLICATE number

  for(irep = 0; irep<((int)NREPS); irep++){

    // Initialize time
    currenttime = 0.0;

    // Initialize the population
      // S on the left-hand side, at carrying capacity K
    for(isite=0; isite< floor(NSITES/2); isite++){
      nW[isite] = (int) 0;
      nS[isite] = (int) floor(K);
    }
      // W on the right-hand side, at K too
      // NB: this is not the carrying capacity of an all-W population when s>0
    for(isite=floor(NSITES/2); isite< NSITES; isite++){
      nW[isite] = (int) floor(K);
      nS[isite] = (int) 0;
    }

    // Compute total population size and total number of infected/uninfected individuals
    totW = 0;
    totS = 0;

    // If only one replicate, print population state
    if(NREPS<2){ // If one sim, print output at regular time intervals
      printPopulationState();
    }

    for(isite=0; isite<NSITES; isite++){
      totS += nS[isite];
      totW += nW[isite];
    }
    totPopSize = totS + totW;

    for(timestep = 0; timestep<((int)MAXSTEPS); timestep++){
      // Check that the population is not extinct, that the drive is still present
      // and that time is not over
      if(totS<1 || totW<1 || currenttime >= MAXTIME){
          break;
        }

      // Do one step of the simulation
      oneStep();

      // If only one replicate, print population state
      if(NREPS<2){ // If one sim, print output at regular time intervals
        if(1.0*timestep - STEPPRINTINTERVAL * floor(timestep/STEPPRINTINTERVAL) < 1.0){
          // if round number of timestep in base STEPPRINTINTERVAL, print ouput
          printPopulationState();
//          printf("Coucou\n");
        }
      }
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
    newnW[isite] = 0;
    newnS[isite] = 0;

    // S individuals
    for(iindiv=0; iindiv<nS[isite]; iindiv++){
      // Compute expected number of offspring
      // per capita (removed * nS[isite])
      fecundity = (1.0 + r * (1.0 - (1.0*(nW[isite] + nS[isite]))/(1.0*K))) * ( 1.0*nS[isite]  / (nW[isite] + nS[isite]) +  1.0*nW[isite] / (nW[isite] + nS[isite]) * wh);
//      printf("B function %.3f\n", (1.0 + r * (1.0 - (1.0*(nW[isite] + nS[isite]))/(1.0*K))));
//      printf("other factor %.3f\n",( 1.0*nS[isite] * nS[isite] / (nW[isite] + nS[isite]) +  nW[isite] * nS[isite] / (nW[isite] + nS[isite]) * wh) );
//      printf("nS = %d", nS[isite]);
//      printf("nW = %d", nW[isite]);
      // Rescale it - can be negative because number of individuals can go over K (migration, independent reproduction)
      if(fecundity < 0.0){
        fecundity = 0.0;
      }
//      printf("FecS = %.3f\n", fecundity);
      deathAndBirth();
      newnS[isite] += newFocal + newOffspring;
    }

    // W individuals
    for(iindiv=0; iindiv<nW[isite]; iindiv++){
      // Compute expected number of offspring in the time interval
      // per capita (removed * nW[isite])
      fecundity = (1.0 + r * (1.0 - (1.0*(nW[isite] + nS[isite]))/(1.0*K))) * ( 1.0*nW[isite]  / (nW[isite] + nS[isite]) * fw + 1.0*nS[isite] / (nW[isite] + nS[isite]) * fw);
      // Rescale it
//      printf("FecW = %.3f\n", fecundity);
      if(fecundity < 0.0){
        fecundity = 0.0;
      }
      deathAndBirth();
      newnW[isite] += newFocal + newOffspring;
    }
  }


  //--------------------------------------------
  // 2) Migration

  // Reset pop
  for(isite=0; isite<NSITES; isite++){
    nW[isite] = 0;
    nS[isite] = 0;
  }

  // Each site, excluding edges
  for(isite=1; isite<(NSITES-1); isite++){
    // W individuals
    nLocal = newnW[isite];
    migrate();
    nW[isite-1] += nLeft;
    nW[isite+1] += nRight;
    nW[isite] += nStay;
    // S individuals
    nLocal = newnS[isite];
    migrate();
    nS[isite-1] += nLeft;
    nS[isite+1] += nRight;
    nS[isite] += nStay;
  }
  // Edges
  // Site 0 -> to the left in fact stay home
    // W individuals
    nLocal = newnW[0];
    migrate();
    nW[0] += nLeft + nStay;
    nW[1] += nRight;
    // S individuals
    nLocal = newnS[0];
    migrate();
    nS[0] += nLeft + nStay;
    nS[1] += nRight;

  // Site (NSITES-1) -> to the right in fact stay home
    // W individuals
    nLocal = newnW[NSITES-1];
    migrate();
    nW[NSITES-1] += nRight + nStay;
    nW[NSITES-2] += nLeft;
    // S individuals
    nLocal = newnS[NSITES-1];
    migrate();
    nS[NSITES-1] += nRight + nStay;
    nS[NSITES-2] += nLeft;

  //--------------------------------------------
  // Updates
  totS = 0;
  totW = 0;

  for(isite=0; isite<NSITES; isite++){
    totW += nW[isite];
    totS += nS[isite];
  }
  totPopSize = totW + totS;
}

/*----------------------------------------------------------------------------*/
/* SUBSTEP: DEATH AND BIRTH */
void deathAndBirth(void){
  // !! Need to have specified fecundity before (depends on genotype)

  // Sample death time
  randnum = drand48(); // Generate uniform random number
  randDeathTime = -1.0/(deathRate) * log(1.0 - randnum); // Convert into exponential random number

  if(randDeathTime < DTIME){ // If the individual dies within DTIME
    deltaTime = randDeathTime; // We are now considering the time until the focal dies
    newFocal = 0; // Individual does not survive
  }else{ // The individual has not died withing DTIME
    deltaTime = DTIME; // Do not change end time of the step
    newFocal = 1; // Individual survives
  };

  newOffspring = 0;
  if(fecundity > 0.0){
    // Sample number of offspring produced during the interval (Poisson)
    // Source: https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
    // (Lambda is below 1 so we are fine with the simple algorithm)
    LPois = exp(-fecundity*deltaTime);
    pPois = 1.0;
    while (pPois >= LPois) {
      newOffspring += 1;
      randnum = drand48();
      pPois *= randnum;
    }
    newOffspring += -1;
  }
  // output is (newFocal, newOffspring)
}

/*----------------------------------------------------------------------------*/
/* SUBSTEP: MIGRATION */
void migrate(void){
  // For a specified site and a specified type, nLocal

  // Reset output variables
  nLeft = 0;
  nRight = 0;
  nStay = 0;

  // Then for each individual... (binomial sampling)
  for(iindiv=0; iindiv<nLocal; iindiv++){
    randnum = drand48(); // Random uniform to decide emigration
    if(randnum < mig){
      // Emigrates
      randnum = drand48(); // Random uniform to pick direction
      if(randnum < 0.5){
        nLeft += 1;
      }else{
        nRight += 1;
      }
    }else{ // Does not emigrate
      nStay += 1;
    }
  }
}

/*----------------------------------------------------------------------------*/
/* IF ONE REP: PRINT POPULATION STATE */
/* PRINT STATE OF THE POPULATION */
void printPopulationState(void){
  // Print current time
  printf("%09.3f, ", currenttime);
  // Print number of individuals of each genotype
  for (isite=0; isite<(NSITES-1); isite++) {
    printf("%04d, ", nS[isite]);
    printf("%04d, ", nW[isite]);
  }
  printf("%04d, ", nS[NSITES-1]);
  printf("%04d\n", nW[NSITES-1]);
}
