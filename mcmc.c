/*
Corona seiirqd covid-19 program
File: mcmc.c (function for mcmc)
Copyright Nicolas Franco - UNamur 2020-2021
Public version 1.0 (28/04/2021)
Version corresponding to the scientific paper:
"Covid-19 Belgium: Extended SEIR-QD model with nursing homes and long-term scenarios-based forecasts"
https://doi.org/10.1101/2020.09.07.20190108
*/

/*
Technical functions for running the mcmc
*/

#include "corona.h"

// create next move with change on all parameters
parameters random_step_all(gsl_rng * r){
    double coef=1;
    parameters par;
  
	par.p0 = gsl_ran_gaussian(r, PO_SGM * coef);
	par.lambdaa = gsl_ran_gaussian(r, LAMBDA_SGM * coef);
	par.lambdas = gsl_ran_gaussian(r, LAMBDA_SGM * coef);
	par.sigma = gsl_ran_gaussian(r, SIGMA_SGM * coef);
	par.tau = gsl_ran_gaussian(r, TAU_SGM * coef);
    par.pay = gsl_ran_gaussian(r, PA_SGM * coef);
    par.pa25 = gsl_ran_gaussian(r, PA_SGM * coef);
    par.pa45 = gsl_ran_gaussian(r, PA_SGM * coef);
    par.pa65 = gsl_ran_gaussian(r, PA_SGM * coef);
    par.pa75 = gsl_ran_gaussian(r, PA_SGM * coef);
    par.pah = gsl_ran_gaussian(r, PA_SGM * coef);
    par.deltay = gsl_ran_gaussian(r, DELTA_SGM * coef);
	par.delta25 = gsl_ran_gaussian(r, DELTA_SGM * coef);
	par.delta45 = gsl_ran_gaussian(r, DELTA_SGM * coef);
	par.delta65 = gsl_ran_gaussian(r, DELTA_SGM * coef);
	par.delta75 = gsl_ran_gaussian(r, DELTA_SGM * coef);
	par.deltah = gsl_ran_gaussian(r, DELTA_SGM * coef);
	par.gammaay = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammaa25 = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammaa45 = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammaa65 = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammaa75 = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammaah = gsl_ran_gaussian(r, GAMMA_SGM * coef);	
	par.gammasy = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammas25 = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammas45 = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammas65 = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammas75 = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammash = gsl_ran_gaussian(r, GAMMA_SGM * coef);
	par.gammaqy = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);
	par.gammaq25 = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);
	par.gammaq45 = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);
	par.gammaq65 = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);
	par.gammaq75 = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);
	par.gammaqh = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);
	par.ry = gsl_ran_gaussian(r, RR_SGM * coef);
	par.r25 = gsl_ran_gaussian(r, RR_SGM * coef);
	par.r45 = gsl_ran_gaussian(r, RR_SGM * coef);
	par.r65 = gsl_ran_gaussian(r, RR_SGM * coef);
	par.r75 = gsl_ran_gaussian(r, RR_SGM * coef);
	par.rh = gsl_ran_gaussian(r, RR_SGM * coef);
	par.rht = gsl_ran_gaussian(r, RR_SGM * coef);
	par.cap = gsl_ran_gaussian(r, CAP_SGM * coef);
	par.caps = gsl_ran_gaussian(r, CAPS_SGM * coef);
	par.pconf = gsl_ran_gaussian(r, PCONF_SGM * coef);
	par.pth = gsl_ran_gaussian(r, PTH_SGM * coef);
	par.pthp = gsl_ran_gaussian(r, PTHP_SGM * coef);
	par.mh = gsl_ran_gaussian(r, MH_SGM * coef);
	par.supphosp = gsl_ran_gaussian(r, PCONF_SGM * coef);
	par.delay = gsl_ran_gaussian(r, DELAY_SGM * coef);
	par.home_lock = gsl_ran_gaussian(r, CONTACT_SGM * coef);
    par.work_lock = gsl_ran_gaussian(r, CONTACT_SGM * coef);
    par.leisure_lock = gsl_ran_gaussian(r, CONTACT_SGM * coef);
	par.home_unlock = gsl_ran_gaussian(r, CONTACT_SGM * coef);
	par.work_unlock = gsl_ran_gaussian(r, CONTACT_SGM * coef);
	par.school_unlock = gsl_ran_gaussian(r, CONTACT_SGM * coef);
	par.recovery = gsl_ran_gaussian(r, RECOVERY_SGM * coef);
	par.rcap = gsl_ran_gaussian(r, DELAY_SGM * coef);
	par.rcaps = gsl_ran_gaussian(r, DELAY_SGM * coef);
	par.leisure_june =  gsl_ran_gaussian(r, CONTACT_SGM * coef);
	par.leisure_july = gsl_ran_gaussian(r, CONTACT_SGM * coef);
	par.leisure_augustus = gsl_ran_gaussian(r, CONTACT_SGM * coef);
	par.school_sept = gsl_ran_gaussian(r, NEWCONTACT_SGM * coef);  
	par.reimp = gsl_ran_gaussian(r, REIMP_SGM * coef);
	par.leisure_sept = gsl_ran_gaussian(r, NEWCONTACT_SGM * coef);
    return par;
}

// create next move with change on one random parameter only

parameters random_step_one(gsl_rng * r){
    double coef=1;
    parameters par;
    int looppar,num;
    
    for (looppar = 1; looppar <= NUMPAR; looppar++) {
		*foreachpar(looppar, &par)=0;		
	}
	
	num=gsl_rng_get(r) % NUMPAR;
	
    switch(num+1){
    case 1:par.p0 = gsl_ran_gaussian(r, PO_SGM * coef);return par;
    case 2:par.lambdaa = gsl_ran_gaussian(r, LAMBDA_SGM * coef);return par;
    case 3:par.lambdas = gsl_ran_gaussian(r, LAMBDA_SGM * coef);return par;
	case 4:par.sigma = gsl_ran_gaussian(r, SIGMA_SGM * coef);return par;
	case 5:par.tau = gsl_ran_gaussian(r, TAU_SGM * coef);return par;
    case 6: par.pay = gsl_ran_gaussian(r, PA_SGM * coef);return par;
    case 7:par.pa25 = gsl_ran_gaussian(r, PA_SGM * coef);return par;
    case 8: par.pa45 = gsl_ran_gaussian(r, PA_SGM * coef);return par;
    case 9: par.pa65 = gsl_ran_gaussian(r, PA_SGM * coef);return par;
    case 10: par.pa75 = gsl_ran_gaussian(r, PA_SGM * coef);return par;
    case 11: par.pah = gsl_ran_gaussian(r, PA_SGM * coef);return par;
    case 12: par.deltay = gsl_ran_gaussian(r, DELTA_SGM * coef);return par;
	case 13:par.delta25 = gsl_ran_gaussian(r, DELTA_SGM * coef);return par;
	case 14:par.delta45 = gsl_ran_gaussian(r, DELTA_SGM * coef);return par;
	case 15:par.delta65 = gsl_ran_gaussian(r, DELTA_SGM * coef);return par;
	case 16:par.delta75 = gsl_ran_gaussian(r, DELTA_SGM * coef);return par;
	case 17:par.deltah = gsl_ran_gaussian(r, DELTA_SGM * coef);return par;
	case 18:par.gammaay = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 19:par.gammaa25 = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 20:par.gammaa45 = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 21:par.gammaa65 = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 22:par.gammaa75 = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 23:par.gammaah = gsl_ran_gaussian(r, GAMMA_SGM * coef);	return par;
	case 24:par.gammasy = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 25:par.gammas25 = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 26:par.gammas45 = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 27:par.gammas65 = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 28:par.gammas75 = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 29:par.gammash = gsl_ran_gaussian(r, GAMMA_SGM * coef);return par;
	case 30:par.gammaqy = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);return par;
	case 31:par.gammaq25 = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);return par;
	case 32:par.gammaq45 = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);return par;
	case 33:par.gammaq65 = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);return par;
	case 34:par.gammaq75 = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);return par;
	case 35:par.gammaqh = gsl_ran_gaussian(r, GAMMAQ_SGM * coef);return par;
	case 36:par.ry = gsl_ran_gaussian(r, RR_SGM * coef);return par;
	case 37:par.r25 = gsl_ran_gaussian(r, RR_SGM * coef);return par;
	case 38:par.r45 = gsl_ran_gaussian(r, RR_SGM * coef);return par;
	case 39:par.r65 = gsl_ran_gaussian(r, RR_SGM * coef);return par;
	case 40:par.r75 = gsl_ran_gaussian(r, RR_SGM * coef);return par;
	case 41:par.rh = gsl_ran_gaussian(r, RR_SGM * coef);return par;
	case 42:par.rht = gsl_ran_gaussian(r, RR_SGM * coef);return par;
	case 43:par.cap = gsl_ran_gaussian(r, CAP_SGM * coef);return par;
	case 44:par.caps = gsl_ran_gaussian(r, CAPS_SGM * coef);return par;
	case 45:par.pconf = gsl_ran_gaussian(r, PCONF_SGM * coef);return par;
	case 46:par.pth = gsl_ran_gaussian(r, PTH_SGM * coef);return par;
	case 47:par.pthp = gsl_ran_gaussian(r, PTHP_SGM * coef);return par;
	case 48:par.mh = gsl_ran_gaussian(r, MH_SGM * coef);return par;
	case 49:par.supphosp = gsl_ran_gaussian(r, PCONF_SGM * coef);return par;
	case 50:par.delay = gsl_ran_gaussian(r, DELAY_SGM * coef);return par;
	case 51:par.home_lock = gsl_ran_gaussian(r, CONTACT_SGM * coef);return par;
	case 52:par.work_lock = gsl_ran_gaussian(r, CONTACT_SGM * coef);return par;
	case 53:par.leisure_lock = gsl_ran_gaussian(r, CONTACT_SGM * coef);return par;
	case 54:par.home_unlock = gsl_ran_gaussian(r, CONTACT_SGM * coef);return par;
	case 55:par.work_unlock = gsl_ran_gaussian(r, CONTACT_SGM * coef);return par;
	case 56:par.school_unlock = gsl_ran_gaussian(r, CONTACT_SGM * coef);return par;
	case 57:par.recovery = gsl_ran_gaussian(r, RECOVERY_SGM * coef);return par;
	case 58:par.rcap = gsl_ran_gaussian(r, DELAY_SGM * coef);return par;
	case 59:par.rcaps = gsl_ran_gaussian(r, DELAY_SGM * coef);return par;
	case 60:par.leisure_june =  gsl_ran_gaussian(r, CONTACT_SGM * coef); return par;
	case 61:par.leisure_july = gsl_ran_gaussian(r, CONTACT_SGM * coef);return par;
	case 62:par.leisure_augustus = gsl_ran_gaussian(r, CONTACT_SGM * coef);return par;
	case 63:par.school_sept = gsl_ran_gaussian(r, NEWCONTACT_SGM * coef);return par;
	case 64:par.reimp = gsl_ran_gaussian(r, REIMP_SGM * coef);return par;
	case 65:par.leisure_sept = gsl_ran_gaussian(r, NEWCONTACT_SGM * coef);return par;
    }
	
    return par;
}

// perform next move (no optimisation)
parameters randomize(parameters temp_par, gsl_rng * r, double coef){
    parameters par;
    if (coef == 0) {
	par = temp_par;
    } else {
    parameters par_step;
    int looppar;
    
    par_step=random_step_all(r);
    for (looppar = 1; looppar <= NUMPAR; looppar++) {
		*foreachpar(looppar, &par)=*foreachpar(looppar, &temp_par)+(*foreachpar(looppar, &par_step) * coef);		
	}

    }
    return par;
}

// perform next move (with optimisation)
parameters opti_randomize(parameters parburninsamples, SEIIRQD * data, SEIIRQD * data_incid, double *burninlikelihood, int burnin, gsl_rng * r,double reimp[4][DAY_MAX + 1]){
    parameters  parbn,par_step,testpar;
	SEIIRQD etatbn[DAY_MAX + 1], etatbn_incid[DAY_MAX + 1];
	int accept,looppar,i;
	double newlikelihoodbn=0,testlikelihoodbn=0,coef;
	
		accept = 0;
		while (accept == 0) {
		    if (burnin == 0) {
		      		  par_step=random_step_all(r);
			coef=COEF_SGM_START_BURNIN;
		    } else if (burnin%2==0){
		      	par_step=random_step_all(r);
		      	coef=1;
		    } else {
		      		  par_step=random_step_one(r);
		      		  coef=COEF_SGM_BURNIN;

		    }
    		for (looppar = 1; looppar <= NUMPAR; looppar++) {
		      *foreachpar(looppar, &parbn)=*foreachpar(looppar, &parburninsamples)+(*foreachpar(looppar, &par_step) * coef);		
	        }
		    accept = test_par(parbn);
		}
		
		newlikelihoodbn = simulikelihood(etatbn, etatbn_incid, data, data_incid,&parbn,r,reimp);		
		
		if (newlikelihoodbn < *burninlikelihood){ 
		*burninlikelihood = newlikelihoodbn;
	
		accept = 1;
		while (accept == 1) {  // search for elarging step
        for (looppar = 1; looppar <= NUMPAR; looppar++) {
		      *foreachpar(looppar, &testpar)=*foreachpar(looppar, &parbn)+(*foreachpar(looppar, &par_step) * coef);		
	     }
	     

		    accept = test_par(testpar);
		   if(accept==1){

		      testlikelihoodbn = simulikelihood(etatbn, etatbn_incid, data, data_incid,&testpar,r,reimp);		

		if (testlikelihoodbn < *burninlikelihood ) {
		*burninlikelihood = testlikelihoodbn;
		  parbn=testpar;
		} else {
		  accept = 0;
		}		      
		   }
		}
		
		// search optimize (step/2)
  
		for(i=0;i<OPTIMIZE;i++){
		  coef = coef/2;
		  // look before
		   for (looppar = 1; looppar <= NUMPAR; looppar++) {
		      *foreachpar(looppar, &testpar)=*foreachpar(looppar, &parbn)-(*foreachpar(looppar, &par_step) * coef);		
	     }
		    accept = test_par(testpar);
		if(accept==1){
		      testlikelihoodbn = simulikelihood(etatbn, etatbn_incid, data, data_incid,&testpar,r,reimp);	
		  if (testlikelihoodbn < *burninlikelihood ) {
		      *burninlikelihood = testlikelihoodbn;
		      parbn=testpar;
		      continue;
		  }
		}
		  // look after
		   for (looppar = 1; looppar <= NUMPAR; looppar++) {
		      *foreachpar(looppar, &testpar)=*foreachpar(looppar, &parbn)+(*foreachpar(looppar, &par_step) * coef);		
	     }
		    accept = test_par(testpar);
		if(accept==1){
		      testlikelihoodbn = simulikelihood(etatbn, etatbn_incid, data, data_incid,&testpar,r,reimp);	
		  if (testlikelihoodbn < *burninlikelihood ) {
		      *burninlikelihood = testlikelihoodbn;
		      parbn=testpar;
		      continue;
		  }
		}		
	  
		}  
		  		   
#if JUST_COMPUTE == 0
		    fprintf(stderr, "              lklh... (opti) %.0f  (%.0f%%)  \r", *burninlikelihood,100.0 * (burnin) / NUM_BURNIN);
#endif
		    return parbn;		    
		} 
    
    // test reverse
    
	coef=-COEF_SGM_BURNIN;


    		for (looppar = 1; looppar <= NUMPAR; looppar++) {
		      *foreachpar(looppar, &parbn)=*foreachpar(looppar, &parburninsamples)+(*foreachpar(looppar, &par_step) * coef);		
	        }
		    accept = test_par(parbn);
		
		
		if(accept==1){
		
		newlikelihoodbn = simulikelihood(etatbn, etatbn_incid, data, data_incid,&parbn,r,reimp);		
    
    if (newlikelihoodbn < *burninlikelihood){ 

		*burninlikelihood = newlikelihoodbn;
		
		accept = 1;
		while (accept == 1) {  // search for elarging step
        for (looppar = 1; looppar <= NUMPAR; looppar++) {
		      *foreachpar(looppar, &testpar)=*foreachpar(looppar, &parbn)+(*foreachpar(looppar, &par_step) * coef);		
	     }
		    accept = test_par(testpar);
		   if(accept==1){
		      testlikelihoodbn = simulikelihood(etatbn, etatbn_incid, data, data_incid,&testpar,r,reimp);		
		if (testlikelihoodbn < *burninlikelihood ) {
		*burninlikelihood = testlikelihoodbn;
		  parbn=testpar;
		} else {
		  accept = 0;
		}		      
		   }
		}
		
		// search optimize (step/2)
		  
		for(i=0;i<OPTIMIZE;i++){
		  coef = coef/2;
		  // look before
		   for (looppar = 1; looppar <= NUMPAR; looppar++) {
		      *foreachpar(looppar, &testpar)=*foreachpar(looppar, &parbn)-(*foreachpar(looppar, &par_step) * coef);		
	     }
		    accept = test_par(testpar);
		if(accept==1){
		      testlikelihoodbn = simulikelihood(etatbn, etatbn_incid, data, data_incid,&testpar,r,reimp);	
		  if (testlikelihoodbn > *burninlikelihood ) {
		      *burninlikelihood = testlikelihoodbn;
		      likelihood(etatbn, etatbn_incid, data, data_incid,&testpar);
		      parbn=testpar;
		      continue;
		  }
		}
		  // look after
		   for (looppar = 1; looppar <= NUMPAR; looppar++) {
		      *foreachpar(looppar, &testpar)=*foreachpar(looppar, &parbn)+(*foreachpar(looppar, &par_step) * coef);		
	     }
		    accept = test_par(testpar);
		if(accept==1){
		      testlikelihoodbn = simulikelihood(etatbn, etatbn_incid, data, data_incid,&testpar,r,reimp);	
		  if (testlikelihoodbn > *burninlikelihood ) {
		      *burninlikelihood = testlikelihoodbn;
		      parbn=testpar;
		      continue;
		  }
		}		
	  
		}  
		  		   
#if JUST_COMPUTE == 0
		    fprintf(stderr, "              lklh... (opti) %.0f  (%.0f%%)  \r", *burninlikelihood,100.0 * (burnin) / NUM_BURNIN);
#endif

		    return parbn;		    
		} 
		}
    
        return parburninsamples;
}


