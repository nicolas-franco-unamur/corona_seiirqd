/*
Corona seiirqd covid-19 program
File: corona_seiirqd.c (main)
Copyright Nicolas Franco - UNamur 2020-2021
Public version 1.0 (28/04/2021)
Version corresponding to the scientific paper:
"Covid-19 Belgium: Extended SEIR-QD model with nursing homes and long-term scenarios-based forecasts"
https://doi.org/10.1101/2020.09.07.20190108
*/

/*
Main file lauching the different steps of the program
*/


#include "corona.h"


// *************************************************************************************

// fct main
int main(int argc, char **argv)
{
    SEIIRQD data[DAY_MAX + 1], data_incid[DAY_MAX + 1];
    double reimp[4][DAY_MAX + 1];
    int looppar;
    const gsl_rng_type *T;
    gsl_rng *r;
    
// initialisation
    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, random_seed());
    readdata(data, data_incid,reimp);

//  Monte-Carlo Metropolis algorithm *************************************************************************************
#if JUST_PLOT_MCMC == 0
    int samplebn;

#if JUST_COMPUTE == 1
    FILE *parameters_mcmc_write;
    printf("%d\n", SIZE_SAMPLE);
#else
    FILE *parameters_mcmc_write;
    parameters_mcmc_write = fopen("inputparameters.txt", "w+");
    printf("%d\n", SIZE_SAMPLE);
    fprintf(parameters_mcmc_write, "%d\n", SIZE_SAMPLE);
#endif

#if READ_PRIORS ==1		// read priors from file
    parameters priors_array[NUM_PRIORS];
    FILE *priors;
    int prior, num_prior;
    priors = fopen("priors.txt", "r");
    fscanf(priors, "%d", &num_prior);
    if (num_prior < NUM_PRIORS) {
	printf("error size priors: %d < %d\n", num_prior, NUM_PRIORS);
	exit(-1);
    }
    for (prior = 0; prior < NUM_PRIORS; prior++) { 
	fscanf(priors, "%d ", &num_prior);
	for (looppar = 1; looppar <= NUMPAR; looppar++) {      
	    fscanf(priors, "%lf", foreachpar(looppar, &priors_array[prior]));
	}
    }
    fclose(priors);
#endif

// start burning phase
    for (samplebn = 0; samplebn < SIZE_SAMPLE_BURNIN; samplebn++) {

	int sample = 0, badburnin = 1, burnin = 0, accept = 0;
	double  burninlikelihood = 0, burninlikelihoodcopy = 0, bnlklh0 = 0;
	parameters parburninsamples;

	while (badburnin) {	// to be sure to be in an admissible zone
	   burninlikelihood = 0;
	   burninlikelihoodcopy = 0;
	   bnlklh0 = 0;
	   
    // if no informed prior:
	    parburninsamples.p0 = PO;
	    parburninsamples.lambdaa = LAMBDA/2.0;
	    parburninsamples.lambdas = LAMBDA;
	    parburninsamples.sigma = SIGMA;
	    parburninsamples.tau = TAU;
	    parburninsamples.pay = PA+2*PA_DELTA;
	    parburninsamples.pa25 = PA+PA_DELTA;
	    parburninsamples.pa45 = PA;
	    parburninsamples.pa65 = PA-PA_DELTA;
	    parburninsamples.pa75 = PA-2*PA_DELTA;
	    parburninsamples.pah = PA-3*PA_DELTA;
	    parburninsamples.deltay = DELTA-2*DELTA_DELTA;
	    parburninsamples.delta25 = DELTA-DELTA_DELTA;
	    parburninsamples.delta45 = DELTA;
	    parburninsamples.delta65 = DELTA+DELTA_DELTA;
	    parburninsamples.delta75 = DELTA+2*DELTA_DELTA;
	    parburninsamples.deltah = DELTA+3*DELTA_DELTA;
	    parburninsamples.gammaay = GAMMA+2*GAMMA_DELTA;
	    parburninsamples.gammaa25 = GAMMA+GAMMA_DELTA;
	    parburninsamples.gammaa45 = GAMMA;
	    parburninsamples.gammaa65 = GAMMA-GAMMA_DELTA;
	    parburninsamples.gammaa75 = GAMMA-2*GAMMA_DELTA;
	    parburninsamples.gammaah = GAMMA-3*GAMMA_DELTA;
	    parburninsamples.gammasy = GAMMA+2*GAMMA_DELTA;
	    parburninsamples.gammas25 = GAMMA+GAMMA_DELTA;
	    parburninsamples.gammas45 = GAMMA;
	    parburninsamples.gammas65 = GAMMA-GAMMA_DELTA;
	    parburninsamples.gammas75 = GAMMA-2*GAMMA_DELTA;
	    parburninsamples.gammash = GAMMA-3*GAMMA_DELTA;
	    parburninsamples.gammaqy = GAMMAQ+2*GAMMAQ_DELTA;
	    parburninsamples.gammaq25 = GAMMAQ+GAMMAQ_DELTA;
	    parburninsamples.gammaq45 = GAMMAQ;
	    parburninsamples.gammaq65 = GAMMAQ-GAMMAQ_DELTA;
	    parburninsamples.gammaq75 = GAMMAQ-2*GAMMAQ_DELTA;
	    parburninsamples.gammaqh = GAMMAQ-3*GAMMAQ_DELTA;
	    parburninsamples.ry = RR-2*RR_DELTA;
	    parburninsamples.r25 = RR-RR_DELTA;
	    parburninsamples.r45 = RR;
	    parburninsamples.r65 = RR+RR_DELTA;
	    parburninsamples.r75 = RR+2*RR_DELTA;
	    parburninsamples.rh = RR+3*RR_DELTA;
	    parburninsamples.rht = RR;
	    parburninsamples.cap = CAP;
	    parburninsamples.caps = CAPS;
	    parburninsamples.pconf = PCONF;
	    parburninsamples.pth = PTH;
	    parburninsamples.pthp = PTHP;
	    parburninsamples.mh = MH;
	    parburninsamples.supphosp = SUPPHOSP;
	    parburninsamples.delay = DELAY;
	    parburninsamples.home_lock = CONTACT_HOME;
	    parburninsamples.work_lock = CONTACT;
	    parburninsamples.leisure_lock = CONTACT;
	    parburninsamples.home_unlock = CONTACT_HOME+CONTACT_DELTA;
	    parburninsamples.work_unlock = CONTACT+CONTACT_DELTA;
	    parburninsamples.school_unlock = CONTACT+CONTACT_DELTA;
	    parburninsamples.recovery = RECOVERY;
	    parburninsamples.rcap = MAX_DAYS_DATA;
	    parburninsamples.rcaps = DELAY;
	    parburninsamples.leisure_june = CONTACT+CONTACT_DELTA;
	    parburninsamples.leisure_july = CONTACT+4*CONTACT_DELTA;
	    parburninsamples.leisure_augustus = CONTACT+2*CONTACT_DELTA;
	    parburninsamples.school_sept = CONTACT+2*CONTACT_DELTA;
	    parburninsamples.reimp = REIMP;
	    parburninsamples.leisure_sept = CONTACT+3*CONTACT_DELTA;

    // if informed prior:
#if READ_PRIORS ==1
	    parburninsamples = priors_array[gsl_rng_get(r) % NUM_PRIORS]; 
#endif

burninlikelihood = DBL_MAX;

for (burnin = 0; burnin < NUM_BURNIN; burnin++) {


// seach with optimized algorithm
#if OPTIMIZE > 0	
    parburninsamples=opti_randomize(parburninsamples, data, data_incid, &burninlikelihood,burnin,r,reimp);
// seach without optimized algorithm (simple random walk)
#else
    parameters parbn;
	SEIIRQD etatbn[DAY_MAX + 1], etatbn_incid[DAY_MAX + 1];
		accept = 0;
		while (accept == 0) {
		    if (burnin == 0) {
			parbn = randomize(parburninsamples, r, COEF_SGM_START_BURNIN);
		    } else {
			parbn = randomize(parburninsamples, r, COEF_SGM_BURNIN);
		    }
		    accept = test_par(parbn);
		}
		newlikelihoodbn = simulikelihood(etatbn, etatbn_incid, data, data_incid,&parbn,r,reimp);				
		if (newlikelihoodbn < burninlikelihood){ 		  

#if JUST_COMPUTE == 0
		    fprintf(stderr, "        lklh... %.0f  (%.0f%%)                \r",  newlikelihoodbn,100.0 * (burnin) / NUM_BURNIN);
#endif

		    burninlikelihood = newlikelihoodbn;
		    parburninsamples = parbn;
		}
		
#endif

	    } // end loop burnin

	    badburnin = 0;
	    if (burninlikelihood > MAX_BURNIN) {
		printf(" restart burn-in      %.2f           \n",  burninlikelihood); 
		badburnin = 1;
	    }
	}
	burninlikelihoodcopy = burninlikelihood;

// end of burnin phase   ****************************************************
#if JUST_COMPUTE == 0
	fprintf(stderr, "\n");
#endif

parameters temp_par;
temp_par = parburninsamples;

	for (sample = 0; sample < SIZE_SAMPLE_PER_BURNIN; sample++) {
	    SEIIRQD etat[DAY_MAX + 1], etat_incid[DAY_MAX + 1];
	    parameters par;
	    double finallikelihood = 0, newlikelihood = 0;
	    int ite = 0;    

	    finallikelihood = simulikelihood(etat, etat_incid, data, data_incid,&temp_par,r,reimp);

	    for (ite = 0; ite < NUM_ITE; ite++) {
		accept = 0;
		while (accept == 0) {
		    par = randomize(temp_par, r, 1);
		    accept = test_par(par);
		}

		newlikelihood = simulikelihood(etat, etat_incid, data, data_incid,&par,r,reimp);

// if best fit FC-HC random walk
#if BEST_FIT == 1
		if (newlikelihood < finallikelihood) {	
// if MCMC
#else

		if (newlikelihood < finallikelihood ||  gsl_rng_uniform(r) < exp(- newlikelihood + finallikelihood)) {	
#endif

		    finallikelihood = newlikelihood;
		    temp_par = par;
		    bnlklh0 = likelihood(etat, etat_incid, data, data_incid,&par);
		    #if JUST_COMPUTE == 0
		    fprintf(stderr, "                lklh sample %.0f    (%.0f%%)               \r",bnlklh0,100.0 * (ite) / NUM_ITE);

#endif
		}

	    }			// end for ite

// print results in file
	    char buffer[2000], temp_buffer[30];
	    sprintf(buffer, "%d ", sample + 1);
	    for (looppar = 1; looppar <= NUMPAR; looppar++) {
		sprintf(temp_buffer, "%.15f ", *foreachpar(looppar, &temp_par));
		strcat(buffer, temp_buffer);
	    }
	    printf("likelihood:%.0f -> %.0f (%.8f)        \n",  burninlikelihoodcopy,bnlklh0, finallikelihood);
	    printf("%s\n", buffer);

#if JUST_COMPUTE == 1
	    parameters_mcmc_write = fopen(OUTPUT_FILE, "a");
#endif

#if BEST_FIT == 1 && JUST_COMPUTE == 1
	    fprintf(parameters_mcmc_write, "likelihood:%.0f -> %.0f\n%s\n", burninlikelihoodcopy,  bnlklh0, buffer);
#else
	    fprintf(parameters_mcmc_write, "%s\n", buffer);
#endif

#if JUST_COMPUTE == 1
	    fclose(parameters_mcmc_write);
#endif
	}			// end for sample

    }				// end looppar burnin

#if JUST_COMPUTE == 1
    return 0;
#else
    fclose(parameters_mcmc_write);
#endif

#endif				//   end just plot 0 = end of MCMC     *****************************************************************

//   output part     *****************************************************************

#if JUST_COMPUTE ==0
output(r,reimp);
#endif				

    return 0;

}
