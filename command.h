/*
Corona seiirqd covid-19 program
File: command.h (with control constants)
Copyright Nicolas Franco - UNamur 2020-2021
Public version 1.0 (28/04/2021)
Version corresponding to the scientific paper:
"Covid-19 Belgium: Extended SEIR-QD model with nursing homes and long-term scenarios-based forecasts"
https://doi.org/10.1101/2020.09.07.20190108
*/

/* *************************************************
This file contains all the elements needed to control the algorithm: to be modified according the desired performances
************************************************* */

#ifndef COMMAND_H_
#define COMMAND_H_

// Define the starting day of the simulation (day1 = 2020/03/01 reported data = 2020/02/29 real situation)
#define DAY_START 1
// Define the end day of the simulation (>=397 accoding to current outputs)
#define DAY_MAX 397
// Define the maximal data day (currently 244)
#define MAX_DAYS_DATA 244

// If define to 1, only perform MCMC with only parameters file inputparameters.txt as output (to be used for parallel programming)
#define JUST_COMPUTE 0    
// If define to 1 (and JUST_COMPUTE=0), just perfom simu from a given parameters file (inputparameters.txt)
#define JUST_PLOT_MCMC 0   
// Name of the output file for parameters  (to be used for parallel programming)
#define OUTPUT_FILE "out.txt"
// to simulate different scenarios (using simucontinu2 instead of simucontinu)
#define PER_CHANGE 0
// if 1, then change and change2 arguments are read as argument to the binary file
#define WITH_ARG 0
// the give a specific name to the outputs (files in plotcs and result_parameters file)
#define PARTNAME "name"
// If several scenarios must be run using the change and change2 parameters (loops)
#define PER_CHANGE_FROM 0
#define PER_CHANGE_TO 0.01
#define PER_CHANGE_STEP 10
#define PER_CHANGE2_FROM 1
#define PER_CHANGE2_TO  1.01
#define PER_CHANGE2_STEP 1

// General speed of the algorithm (1 => 2000 nursing home, >1 => approximation with less precision)
#define SPEED   2
// Set to >1 if several runs of the stochastic part must be done to get the likelihood (default 1)
#define MOYRAND 1

// If informed priors must be read from priors.txt file
#define READ_PRIORS 1
// num of those priors
#define NUM_PRIORS 2900

// Number of iteration for the bunrning phase (at least 20000000 for initial calibration  but this can take days / >=2000 from informed priors)
#define NUM_BURNIN 2000	
// Number of runs of bunning (= restart)
#define SIZE_SAMPLE_BURNIN 1
// Maximal value of the likelihood (if not reached, then restart a burning)
#define MAX_BURNIN 0
// >1 if optimisation seach must be used (default 5)
#define OPTIMIZE 5  

// Number of iteration between each sample (default 20000)
#define NUM_ITE 200		
// Number of samples per burning (default 10)
#define SIZE_SAMPLE_PER_BURNIN 10	
// If best-fit mode must be used instead of MCMC-Hasting (default 0)
#define BEST_FIT 0

// Step multiplicative coefficient for the burning perriod (typically 5)
#define COEF_SGM_BURNIN 5	
// Coefficient of the normal distribution at the start of the burning period (tipically up to 100)
#define COEF_SGM_START_BURNIN 100	
// Coefficient of the normal distribution at the start of the mcmc  (default 0)
#define COEF_SGM_START 0	
// not to change
#define SIZE_SAMPLE (SIZE_SAMPLE_BURNIN*SIZE_SAMPLE_PER_BURNIN)
// malus to the likelihood in cas of no admissible point
#define MALUS_LIKELIHOOD 10000000

#endif