// To compile:
// gcc stochsim_base.c -lm -o test

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h> // to get system time

// Simulation parameters
#define DTIME 0.05 // Time interval at which population is updated (has to be <1)

#define MAXSTEPS 10000000000 // Maximum number of simulation steps (in case a loop goes wrong)
#define STEPPRINTINTERVAL 100 // Print result every ... steps


// Define variables
int nD; // number of DD individuals
int nO; // number of OO individuals
// Note: WT genotype written with letter O
int newnD; // temporary number of DD
int newnO; // temporary number of OO
int newFocal; // Whether the individual that we are considering survives
int newOffspring; // Number of offspring produced by the focal
double fecundityD; // Expected number of offspring in the time interval
double fecundityO; // Expected number of offspring in the time interval
double fecundity;

int timestep; // Time Step
double currenttime; // Time value
double deltaTime; // Time interval, if death occurs within DTIME

int iindiv; // Counter
double randnum; // Random number
double randDeathTime;
double randBirthTime;

double pPois; // Temp variable for Poisson sampling
double LPois; // Temp variable for Poisson sampling

// Functions used in the script
void printPopulationState(void); // Print current state of the population
void oneStep(void); // Do one step of the simulation and update the population accordingly
void deathAndBirth(void); // substep: Death and Birth

int g_argc; // Global variable for command line parameter entry
char **g_argv; // Global variable for command line parameter entry

// Introduce parameter names (values entered when running the script)
double K; // Carrying capacity
double s; // Zygote survival cost induced by the drive (0<= s <= 1)
double r; // Growth rate
int N0D; // Introduced number of Drive individuals
double MAXTIME; // Maximum time until which simulation is run
int NREPS; // Number of simulation replicates

// Other parameters
double deathRate = 1.0; // Scaled parameter

/* Function to process the parameters, give them names */
void process_command_line(void)
{
  // Extract parameters from the command line
  /// Model parameters
  K = (double) strtod(g_argv[1], (char **) NULL); // K, carrying capacity
  s = (double) strtod(g_argv[2], (char **) NULL); // Zygote survival cost induced by the drive (0<= s <=1)
  r = (double) strtod(g_argv[3], (char **) NULL); // Growth rate
  N0D = (int) strtod(g_argv[4], (char **) NULL); // Number of introduced D individuals
  /// Simulation parameters
  MAXTIME = (double) strtod(g_argv[5], (char **) NULL); // Maximal time value of the simulation
  NREPS = (int) strtod(g_argv[6], (char **) NULL); // Number of replicates
}

int main(int argc, char **argv)
{
  // Process command line arguments
  g_argc = argc;
  g_argv = argv;
  process_command_line();
  printf("K = %.0f, s = %.4f, r = %.2f, N0D = %d, MAXTIME = %.2f, NREPS = %d\n", K, s, r, N0D, MAXTIME, NREPS);

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
    nO = K - N0D;
    nD = N0D;
    
    // If only one replicate, print population state
    if(NREPS<2){ // If one sim, print output at regular time intervals
      printPopulationState();
    }
        
    for(timestep = 0; timestep<((int)MAXSTEPS); timestep++){
      // Check that the population is not extinct, that the drive is still present
      // and that time is not over
      if((nO + nD)<1 || currenttime >= MAXTIME){
          break;
        }
      
      // Do one step of the simulation
      oneStep();

      // If only one replicate, print population state
      if(NREPS<2){ // If one sim, print output at regular time intervals
        if(1.0*timestep - STEPPRINTINTERVAL * floor(timestep/STEPPRINTINTERVAL) < 1.0){
          // if round number of timestep in base STEPPRINTINTERVAL, print ouput
          printPopulationState();
        }
      }
    }
    
    // Simulation finished -> save final outcome
    if(NREPS>1){
    // Save final population step
      // TODO THINK ABOUT WHAT TO PRINT WHEN MULTIPLE REPLICATES
    }else{
      printPopulationState();
    }
  } //end rep
  
  return(0);

} // end main Functions



/*----------------------------------------------------------------------------*/
/* ONE SIMULATION STEP - ZYGOTE CONVERSION */
void oneStep(void){
  // Update time
  currenttime += (1.0)*DTIME;
//  printf("\n");
  
  //-------------------------------------------
  // 1) Death and birth
    // Reset newpop
    newnD = 0;
    newnO = 0;

  int totNewOffO = 0;
    // WT individuals
    for(iindiv=0; iindiv<nO; iindiv++){
//      printf("coucou\n");
      // Compute expected number of offspring in the time interval
      fecundityO = (1.0 + r * (1.0 - (1.0*(nD + nO))/(1.0*K))) * (nO) / (nO + nD) * deltaTime;
      fecundity = fecundityO;
      // Rescale it - can be negative because number of individuals can go over K (migration, independent reproduction)
      if(fecundity < 0.0){
        fecundity = 0.0;
      }
      deathAndBirth();
      newnO += newFocal + newOffspring;
      totNewOffO += newOffspring;
    }
    
  int totNewOffD = 0;
    // Drive individuals
    for(iindiv=0; iindiv<nD; iindiv++){
      // Compute expected number of offspring in the time interval
      fecundityD = (1.0 + r * (1.0 - (1.0*(nD + nO))/(1.0*K))) * (1.0 - s) * (2.0 * nO + nD) / (nO + nD) * deltaTime; // Expected number of offspring in the time interval
      fecundity = fecundityD;
      // Rescale it
      if(fecundity < 0.0){
        fecundity = 0.0;
      }
      deathAndBirth();
      newnD += newFocal + newOffspring;
      totNewOffD += newOffspring;
    }
//  printf("%d, %d\n", totNewOffO, totNewOffD);
  
  
//  //#######################
//  // Print current time
//  printf("%.3f, ", currenttime);
//  // Print number of individuals of each genotype
//  for (isite=0; isite<(NSITES-1); isite++) {
//    printf("%d, ", newnO[isite]);
//    printf("%d, ", newnD[isite]);
//  }
//  printf("%d, ", newnO[NSITES-1]);
//  printf("%d\n", newnD[NSITES-1]);
//  //###########################
  
  


  //--------------------------------------------
  // Updates
  nD = newnD;
  nO = newnO;
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
    newFocal = 1;
  };
  
  newOffspring = 0;
  if(fecundity > 0.0){
    // Sample number of offspring produced during the interval (Poisson)
    // Source: https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
    // (Lambda is below 1 so we are fine with the simple algorithm)
    LPois = exp(-fecundity);
    pPois = 1.0;
    while (pPois >= LPois) {
      newOffspring += 1;
      randnum = drand48();
      pPois *= randnum;
    }
    newOffspring += -1;
  }
//  printf("%d, ", newOffspring);
  
  // output is (newFocal, newOffspring)
//  printf("%0.4f, %.2f, %d, %d\n", (1.0 - (1.0*(nD[isite] + nO[isite]))/(1.0*K)), fecundity, newFocal, newOffspring);
}

/*----------------------------------------------------------------------------*/
/* IF ONE REP: PRINT POPULATION STATE */
/* PRINT STATE OF THE POPULATION */
void printPopulationState(void){
  // Print current time
  printf("%07.3f", currenttime);
  printf(", randDeathTime = %.3f", randDeathTime);
  printf(", fO = %.5f, fD = %.5f, ", fecundityO, fecundityD);
  printf("nO = %04d, nD = %04d\n", nO, nD);
}
