/*
Corona seiirqd covid-19 program
File: output.c (function for generating output)
Copyright Nicolas Franco - UNamur 2020-2021
Public version 1.0 (28/04/2021)
Version corresponding to the scientific paper:
"Covid-19 Belgium: Extended SEIR-QD model with nursing homes and long-term scenarios-based forecasts"
https://doi.org/10.1101/2020.09.07.20190108
*/

/*
Fonction creating outputs of simulation: to be modified with the desired outputs
*/

#include "corona.h"


void output(gsl_rng * r,double reimp[4][DAY_MAX + 1]){
    int looppar;

// read parameters
    parameters parsamples[SIZE_SAMPLE];
    FILE *parameters_mcmc;
    int sample, num_sample, day;
    


    parameters_mcmc = fopen("inputparameters.txt", "r");
    fscanf(parameters_mcmc, "%d", &num_sample);
    if (num_sample < SIZE_SAMPLE) {
	printf("error size: %d < %d\n", num_sample, SIZE_SAMPLE);
	exit(-1);
    }
    for (sample = 0; sample < SIZE_SAMPLE; sample++) {
	fprintf(stderr, " read... %.0f%%                  \r", 100.0 * (sample) / SIZE_SAMPLE);
	fscanf(parameters_mcmc, "%d ", &num_sample);
	for (looppar = 1; looppar <= NUMPAR; looppar++) {
	    fscanf(parameters_mcmc, "%lf", foreachpar(looppar, &parsamples[sample]));
	}
    }
    fclose(parameters_mcmc);
    
// gather results (perform simu until DAY_MAX days)
    int loop;
    double *tempsetatarray = NULL, change = 0,change2=0;
    SEIIRQDarray etatarray[DAY_MAX + 1], etatarray_incid[DAY_MAX + 1];
    SEIIRQD etatsamples[SIZE_SAMPLE][DAY_MAX + 1], etatsamples_incid[SIZE_SAMPLE][DAY_MAX + 1];

#if PER_CHANGE > 0     
    for (change = PER_CHANGE_FROM; change <= PER_CHANGE_TO; change += PER_CHANGE_STEP) {
    for (change2 = PER_CHANGE2_FROM; change2 <= PER_CHANGE2_TO; change2 += PER_CHANGE2_STEP) {
#endif


#if WITH_ARG > 0
if(argc<3){
 printf("error argc\n");
 exit(-1);   
}
change = atof(argv[1]);
change2 = atof(argv[2]);
#endif

	for (sample = 0; sample < SIZE_SAMPLE; sample++) {
	    fprintf(stderr, " launch simu... %.0f%%  (%2.2f-%s)                \r", 100.0 * (sample) / SIZE_SAMPLE, change,PARTNAME);
	    init(etatsamples[sample], etatsamples_incid[sample], &parsamples[sample]);
	    if (change2 > 0) {
		simucontinu2(etatsamples[sample], etatsamples_incid[sample], &(parsamples[sample]), DAY_MAX, change, change2, r,reimp);
	    } else {
		simucontinu(etatsamples[sample], etatsamples_incid[sample], &(parsamples[sample]), DAY_MAX, r,reimp);
	    }

	    for (day = 0; day <= DAY_MAX; day++) {	//transpose
		for (loop = 1; loop <= NUMSEG; loop++) {
		    tempsetatarray = foreachetatarray(loop, &etatarray[day]);
		    tempsetatarray[sample] = *foreachetat(loop, &etatsamples[sample][day]);
		    tempsetatarray = foreachetatarray(loop, &etatarray_incid[day]);
		    tempsetatarray[sample] = *foreachetat(loop, &etatsamples_incid[sample][day]);
		}
	    }

	}			// fin for sample (2)


// new plot (all P5%) *******************************************************
	FILE *result_P5 = NULL;
	int perc;
	double *tempsarray = NULL, percpouc;
	char namefile[50], namefilepar[50];
	parametersarray pararray;
	double rzero[SIZE_SAMPLE];
	double rzerop[SIZE_SAMPLE];
	double rzeropp[SIZE_SAMPLE];
	double rzerotb[SIZE_SAMPLE];
	double rzerotbp[SIZE_SAMPLE];
	double rzerotbpp[SIZE_SAMPLE];
	double rzerote[SIZE_SAMPLE];
	double rzerotep[SIZE_SAMPLE];
	double rzerotepp[SIZE_SAMPLE];

	double rzeromay[SIZE_SAMPLE];
	double rzerojune[SIZE_SAMPLE];
	double rzerojuly[SIZE_SAMPLE];
	double rzeroaug[SIZE_SAMPLE];
	double rzerosept[SIZE_SAMPLE];
	double rzerotmay[SIZE_SAMPLE];
	double rzerotjune[SIZE_SAMPLE];
	double rzerotjuly[SIZE_SAMPLE];
	double rzerotaug[SIZE_SAMPLE];
	double rzerotsept[SIZE_SAMPLE];	
	
		double rzeroh[SIZE_SAMPLE];
		double rzeroht[SIZE_SAMPLE];

	
	for (sample = 0; sample < SIZE_SAMPLE; sample++) {

	    for (looppar = 1; looppar <= NUMPAR; looppar++) {
		tempsarray = foreachpararray(looppar, &pararray);
		tempsarray[sample] = *foreachpar(looppar, &parsamples[sample]);
	    }

	   	rzero[sample] =compute_rzero(&parsamples[sample],1,1,1,1);
	    rzerop[sample] =compute_rzero(&parsamples[sample],1,0,parsamples[sample].leisure_lock,1);
	    rzeropp[sample] =compute_rzero(&parsamples[sample],parsamples[sample].work_lock,0,parsamples[sample].leisure_lock,parsamples[sample].home_lock);
	    
	    rzeromay[sample] =compute_rzero(&parsamples[sample],parsamples[sample].work_unlock,0.6*parsamples[sample].school_unlock,parsamples[sample].leisure_lock,parsamples[sample].home_lock+(parsamples[sample].home_unlock-parsamples[sample].home_lock)/2);
	    rzerojune[sample] =compute_rzero(&parsamples[sample],parsamples[sample].work_unlock,1*parsamples[sample].school_unlock,parsamples[sample].leisure_june,parsamples[sample].home_unlock);
        rzerojuly[sample] =compute_rzero(&parsamples[sample],parsamples[sample].work_unlock,0,parsamples[sample].leisure_july,parsamples[sample].home_unlock);
        rzeroaug[sample] = compute_rzero(&parsamples[sample],parsamples[sample].work_unlock,0,parsamples[sample].leisure_augustus,parsamples[sample].home_unlock);
        rzerosept[sample] =compute_rzero(&parsamples[sample],parsamples[sample].work_unlock,parsamples[sample].school_sept,parsamples[sample].leisure_sept,parsamples[sample].home_unlock);	        

	    rzerotb[sample] = rzero[sample] * etatarray[1].S[sample] / (NPOP - etatarray[1].D[sample]);
	    rzerote[sample] = rzero[sample] * etatarray[14].S[sample] / (NPOP - etatarray[14].D[sample]);
	    rzerotbp[sample] = rzerop[sample] * etatarray[15].S[sample] / (NPOP - etatarray[15].D[sample]);
	    rzerotep[sample] = rzerop[sample] * etatarray[19].S[sample] / (NPOP - etatarray[19].D[sample]);
	    rzerotbpp[sample] = rzeropp[sample] * etatarray[20].S[sample] / (NPOP - etatarray[20].D[sample]);
	    rzerotepp[sample] = rzeropp[sample] * etatarray[65].S[sample] / (NPOP - etatarray[65].D[sample]);
	    rzerotmay[sample] = rzeromay[sample] * etatarray[95].S[sample] / (NPOP - etatarray[95].D[sample]);
	    rzerotjune[sample] = rzerojune[sample] * etatarray[101].S[sample] / (NPOP - etatarray[101].D[sample]);
	    rzerotjuly[sample] = rzerojuly[sample] * etatarray[124].S[sample] / (NPOP - etatarray[124].D[sample]);
	    rzerotaug[sample] = rzeroaug[sample] * etatarray[124].S[sample] / (NPOP - etatarray[124].D[sample]);
	    rzerotsept[sample] = rzerosept[sample] * etatarray[MAX_DAYS_DATA].S[sample] / (NPOP - etatarray[MAX_DAYS_DATA].D[sample]);

	   	rzeroh[sample] =parsamples[sample].mh*parsamples[sample].lambdaa*(parsamples[sample].pah/parsamples[sample].gammaah + (1-parsamples[sample].pah)/parsamples[sample].tau)	+ parsamples[sample].mh*parsamples[sample].lambdas*(1-parsamples[sample].pah)/(parsamples[sample].gammash+parsamples[sample].deltah+ (parsamples[sample].rht*parsamples[sample].pconf)/2);

        rzeroht[sample] = rzeroh[sample] * etatarray[MAX_DAYS_DATA].SH[sample] / (NPOPH);

	        
	}

	for (day = DAY_START; day <= DAY_MAX; day++) {
	    fprintf(stderr, " sort... %.0f%%                  \r", 100.0 * (day) / DAY_MAX);
	    for (loop = 1; loop <= NUMSEG; loop++) {
		tempsetatarray = foreachetatarray(loop, &etatarray[day]);
		gsl_sort(tempsetatarray, 1, SIZE_SAMPLE);
		tempsetatarray = foreachetatarray(loop, &etatarray_incid[day]);
		gsl_sort(tempsetatarray, 1, SIZE_SAMPLE);
	    }
	}

	for (looppar = 1; looppar <= NUMPAR; looppar++) {
	    tempsarray = foreachpararray(looppar, &pararray);
	    gsl_sort(tempsarray, 1, SIZE_SAMPLE);
	}
	gsl_sort(rzero, 1, SIZE_SAMPLE);
	gsl_sort(rzerop, 1, SIZE_SAMPLE);
	gsl_sort(rzeropp, 1, SIZE_SAMPLE);
	gsl_sort(rzerotb, 1, SIZE_SAMPLE);
	gsl_sort(rzerotbp, 1, SIZE_SAMPLE);
	gsl_sort(rzerotbpp, 1, SIZE_SAMPLE);
	gsl_sort(rzerote, 1, SIZE_SAMPLE);
	gsl_sort(rzerotep, 1, SIZE_SAMPLE);
	gsl_sort(rzerotepp, 1, SIZE_SAMPLE);

	gsl_sort(rzeromay, 1, SIZE_SAMPLE);
	gsl_sort(rzerojune, 1, SIZE_SAMPLE);
	gsl_sort(rzerojuly, 1, SIZE_SAMPLE);
	gsl_sort(rzeroaug, 1, SIZE_SAMPLE);
	gsl_sort(rzerosept, 1, SIZE_SAMPLE);
	gsl_sort(rzerotmay, 1, SIZE_SAMPLE);
	gsl_sort(rzerotjune, 1, SIZE_SAMPLE);
	gsl_sort(rzerotjuly, 1, SIZE_SAMPLE);
	gsl_sort(rzerotaug, 1, SIZE_SAMPLE);
	gsl_sort(rzerotsept, 1, SIZE_SAMPLE);
	
		gsl_sort(rzeroh, 1, SIZE_SAMPLE);
		gsl_sort(rzeroht, 1, SIZE_SAMPLE);


	
	for (perc = 5; perc < 100; perc += 5) {
	    percpouc = 0.01 * (double) perc;
	    fprintf(stderr, " percentiles... %.0f%%                  \r", percpouc);
	    result_P5 = NULL;

	    if (change2 > 0) {
#if PARTPLOT ==1
		sprintf(namefile, "plotcs/p%d-%s-c%2.1f.txt", perc,PARTNAME,change);
#else
		sprintf(namefile, "plotcs/p%d-c%2.1f-%1.1f.txt", perc, change,change2);
#endif
	    } else {
		sprintf(namefile, "plot/p%d.txt", perc);

	    }
	    result_P5 = fopen(namefile, "w+");

	    for (day = DAY_START; day <= DAY_MAX; day++) {

		fprintf(result_P5, "%d ", day);
		for (loop = 1; loop <= NUMSEG; loop++) {
		    tempsetatarray = foreachetatarray(loop, &etatarray[day]);
		    fprintf(result_P5, "%.2f ", gsl_stats_quantile_from_sorted_data(tempsetatarray, 1, SIZE_SAMPLE, percpouc));
		}
		for (loop = 1; loop <= NUMSEG; loop++) {
		    tempsetatarray = foreachetatarray(loop, &etatarray_incid[day]);
		    fprintf(result_P5, "%.2f ", gsl_stats_quantile_from_sorted_data(tempsetatarray, 1, SIZE_SAMPLE, percpouc));
		}
		fprintf(result_P5, "\n");
	    }

	    fclose(result_P5);
	}


// plot of specific stats in file result_parameters    ****************************************************
	FILE *result_par = NULL;
	if (change2 > 0) {
#if PARTPLOT ==1
	    sprintf(namefilepar, "result_parameters_%s_cs%2.1f.txt",PARTNAME,change);
#else
	    sprintf(namefilepar, "result_parameters_cs%2.1f-%1.1f.txt", change,change2);
#endif
	} else {
	    sprintf(namefilepar, "result_parameters.txt");
	}
	result_par = fopen(namefilepar, "w+");

	fprintf(result_par, "results on day: %d\n", MAX_DAYS_DATA);
	tempsarray = foreachpararray(1, &pararray);
	fprintf(result_par, "CI: mean=%.12f med=%.12f [%.12f ; %.12f]  \n",
		gsl_stats_mean(tempsarray, 1, SIZE_SAMPLE) * (NPOP - NPOPH),
		gsl_stats_median_from_sorted_data(tempsarray, 1, SIZE_SAMPLE) * (NPOP - NPOPH),
		gsl_stats_quantile_from_sorted_data(tempsarray, 1, SIZE_SAMPLE, 0.05) * (NPOP - NPOPH), gsl_stats_quantile_from_sorted_data(tempsarray, 1, SIZE_SAMPLE, 0.95) * (NPOP - NPOPH));

	for (looppar = 1; looppar <= NUMPAR; looppar++) {
	    tempsarray = foreachpararray(looppar, &pararray);
	    fprintf(result_par, "%s: mean=%.12f med=%.12f [%.12f ; %.12f]  \n",
		    foreachparname(looppar),
		    gsl_stats_mean(tempsarray, 1, SIZE_SAMPLE),
		    gsl_stats_median_from_sorted_data(tempsarray, 1, SIZE_SAMPLE),
		    gsl_stats_quantile_from_sorted_data(tempsarray, 1, SIZE_SAMPLE, 0.05), gsl_stats_quantile_from_sorted_data(tempsarray, 1, SIZE_SAMPLE, 0.95));
	}

fprintf(result_par,"rzero : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzero , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzero , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzero, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzero , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzerop : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerop , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerop , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerop, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerop , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzeropp : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzeropp , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzeropp , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzeropp, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzeropp , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzerotb : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotb , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotb , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotb, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotb , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzerote : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerote , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerote , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerote, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerote , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzerotbp : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotbp , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotbp , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotbp, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotbp , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzerotep : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotep , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotep , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotep, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotep , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzerotbpp : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotbpp , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotbpp , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotbpp, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotbpp , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzerotepp : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotepp , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotepp , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotepp, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotepp , 1, SIZE_SAMPLE,0.95)); 



fprintf(result_par,"rzeromay : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzeromay , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzeromay , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzeromay, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzeromay , 1, SIZE_SAMPLE,0.95));


fprintf(result_par,"rzerotmay : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotmay , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotmay , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotmay, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotmay , 1, SIZE_SAMPLE,0.95));


fprintf(result_par,"rzerojune : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerojune , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerojune , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerojune, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerojune , 1, SIZE_SAMPLE,0.95));


fprintf(result_par,"rzerotjune : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotjune , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotjune , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotjune, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotjune , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzerojuly : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerojuly , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerojuly , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerojuly, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerojuly , 1, SIZE_SAMPLE,0.95));


fprintf(result_par,"rzerotjuly : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotjuly , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotjuly , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotjuly, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotjuly , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzeroaug : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzeroaug , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzeroaug , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzeroaug, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzeroaug , 1, SIZE_SAMPLE,0.95));


fprintf(result_par,"rzerotaug : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotaug , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotaug , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotaug, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotaug , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzerosept : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerosept , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerosept , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerosept, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerosept , 1, SIZE_SAMPLE,0.95));


fprintf(result_par,"rzerotsept : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotsept , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotsept , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotsept, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotsept , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzeroh : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzeroh , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzeroh , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzeroh, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzeroh , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"rzeroht : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzeroht , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzeroht , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzeroht, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzeroht , 1, SIZE_SAMPLE,0.95));

printf("rzerotjune : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotjune , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotjune , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotjune, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotjune , 1, SIZE_SAMPLE,0.95));

printf("rzerojuly : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerojuly , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerojuly , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerojuly, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerojuly , 1, SIZE_SAMPLE,0.95));


printf("rzerotjuly : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotjuly , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotjuly , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotjuly, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotjuly , 1, SIZE_SAMPLE,0.95));

printf("rzeroaug : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzeroaug , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzeroaug , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzeroaug, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzeroaug , 1, SIZE_SAMPLE,0.95));


printf("rzerotaug : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotaug , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotaug , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotaug, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotaug , 1, SIZE_SAMPLE,0.95));

printf("rzerosept : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerosept , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerosept , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerosept, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerosept , 1, SIZE_SAMPLE,0.95));


printf("rzerotsept : mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotsept , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotsept , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotsept, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotsept , 1, SIZE_SAMPLE,0.95));


// additional stats
double stat[SIZE_SAMPLE];

	fprintf(result_par, "\n");

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= etatarray[DAY_MAX].D[sample];
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Deaths on end : mean=%.0f med=%.0f [%.0f ; %.0f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*etatarray[DAY_MAX].D[sample]/(NPOP);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Deaths on end : mean=%.3f%% med=%.3f%% [%.3f%% ; %.3f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= etatarray[MAX_DAYS_DATA].E[sample]+etatarray[MAX_DAYS_DATA].AI[sample]+etatarray[MAX_DAYS_DATA].PI[sample]+etatarray[MAX_DAYS_DATA].SI[sample]+etatarray[MAX_DAYS_DATA].Q[sample];
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Infected on day : mean=%.0f med=%.0f [%.0f ; %.0f]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[MAX_DAYS_DATA].E[sample]+etatarray[MAX_DAYS_DATA].AI[sample]+etatarray[MAX_DAYS_DATA].PI[sample]+etatarray[MAX_DAYS_DATA].SI[sample]+etatarray[MAX_DAYS_DATA].Q[sample])/(NPOP-etatarray[MAX_DAYS_DATA].D[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Infected on day : mean=%.3f%% med=%.3f%% [%.3f%% ; %.3f%%]    \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[MAX_DAYS_DATA].E[sample]+etatarray[MAX_DAYS_DATA].AI[sample]+etatarray[MAX_DAYS_DATA].PI[sample])/(etatarray[MAX_DAYS_DATA].E[sample]+etatarray[MAX_DAYS_DATA].AI[sample]+etatarray[MAX_DAYS_DATA].PI[sample]+etatarray[MAX_DAYS_DATA].SI[sample]+etatarray[MAX_DAYS_DATA].Q[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc Asymptomatic (on day) :  mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= etatarray[MAX_DAYS_DATA].SI[sample];
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"I symptomatic on day : mean=%.0f med=%.0f [%.0f ; %.0f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[MAX_DAYS_DATA].R[sample]-etatarray[MAX_DAYS_DATA].RH[sample])/(NPOP-etatarray[MAX_DAYS_DATA].D[sample]-etatarray[MAX_DAYS_DATA].SH[sample]-etatarray[MAX_DAYS_DATA].EH[sample]-etatarray[MAX_DAYS_DATA].AIH[sample]-etatarray[MAX_DAYS_DATA].PIH[sample]-etatarray[MAX_DAYS_DATA].SIH[sample]-etatarray[MAX_DAYS_DATA].RH[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune (not considering homes) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[MAX_DAYS_DATA].EH[sample]+etatarray[MAX_DAYS_DATA].AIH[sample]+etatarray[MAX_DAYS_DATA].PIH[sample]+etatarray[MAX_DAYS_DATA].SIH[sample])/(etatarray[MAX_DAYS_DATA].SH[sample]+etatarray[MAX_DAYS_DATA].EH[sample]+etatarray[MAX_DAYS_DATA].AIH[sample]+etatarray[MAX_DAYS_DATA].PIH[sample]+etatarray[MAX_DAYS_DATA].SIH[sample]+etatarray[MAX_DAYS_DATA].RH[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Infected on day in homes : mean=%.3f%% med=%.3f%% [%.3f%% ; %.3f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[MAX_DAYS_DATA].RH[sample])/(NPOPH);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune in homes (vs total pop homes) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));




for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma;
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Latent (pre-infectious) period : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/(parsamples[sample].tau);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Pre-symptomatic period  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].tau);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Incubation period for symptomatic (latent + pre-symptomatic)  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].gammaay);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for asymptomatic  0-24  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].gammaa25);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for asymptomatic  25-44  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].gammaa45);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for asymptomatic 45-64  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].gammaa65);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for asymptomatic  65-74  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].gammaa75);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for asymptomatic  75+  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].gammaah);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for asymptomatic  home  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].tau)+1/(parsamples[sample].gammasy);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for symptomatic  0-24  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].tau)+1/(parsamples[sample].gammas25);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for symptomatic  25-44  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].tau)+1/(parsamples[sample].gammas45);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for symptomatic  45-64  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].tau)+1/(parsamples[sample].gammas65);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for symptomatic  65-74  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].tau)+1/(parsamples[sample].gammas75);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for symptomatic  75+  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/parsamples[sample].sigma+1/(parsamples[sample].tau)+1/(parsamples[sample].gammash);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Total disease duration for symptomatic  home  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));





for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/(parsamples[sample].ry+ parsamples[sample].gammaqy);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"hospitalisation duration 0-24  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/(parsamples[sample].r25+ parsamples[sample].gammaq45);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"hospitalisation duration 25-44  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/(parsamples[sample].r45+ parsamples[sample].gammaq45);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"hospitalisation duration 45-64  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/(parsamples[sample].r65+ parsamples[sample].gammaq65);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"hospitalisation duration 65-74  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/(parsamples[sample].r75+ parsamples[sample].gammaq75);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"hospitalisation duration 75+  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 1/(parsamples[sample].rh+ parsamples[sample].gammaqh);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"hospitalisation duration homes  : mean=%.1f days med=%.1f days [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D[sample])/(etatarray[214].D[sample]+etatarray[214].R[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].DY[sample])/(etatarray[214].DY[sample]+etatarray[214].RY[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 0-24: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D25[sample])/(etatarray[214].D25[sample]+etatarray[214].R25[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 25-44: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D45[sample])/(etatarray[214].D45[sample]+etatarray[214].R45[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 45-64: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D65[sample])/(etatarray[214].D65[sample]+etatarray[214].R65[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 65-74: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D75H[sample])/(etatarray[214].D75H[sample]+etatarray[214].R75[sample]+etatarray[MAX_DAYS_DATA].RH[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 75+ (nurs. homes included): mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[93].D[sample])/(etatarray[93].D[sample]+etatarray[93].R[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR March-April: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[93].DY[sample])/(etatarray[93].DY[sample]+etatarray[93].RY[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 0-24 March-April: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[93].D25[sample])/(etatarray[93].D25[sample]+etatarray[93].R25[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 25-44 March-April: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[93].D45[sample])/(etatarray[93].D45[sample]+etatarray[93].R45[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 45-64 March-April: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[93].D65[sample])/(etatarray[93].D65[sample]+etatarray[93].R65[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 65-74 March-April: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[93].D75H[sample])/(etatarray[93].D75H[sample]+etatarray[93].R75[sample]+etatarray[93].RH[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR 75+ March-April (nurs. homes included): mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D[sample]-etatarray[122].D[sample])/(etatarray[214].D[sample]+etatarray[214].R[sample]-etatarray[122].D[sample]-etatarray[122].R[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR July-September: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].DY[sample]-etatarray[122].DY[sample])/(etatarray[214].DY[sample]+etatarray[214].RY[sample]-etatarray[122].DY[sample]-etatarray[122].RY[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR July-September 0-24: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D25[sample]-etatarray[122].D25[sample])/(etatarray[214].D25[sample]+etatarray[214].R25[sample]-etatarray[122].D25[sample]-etatarray[122].R25[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR July-September 25-44: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D45[sample]-etatarray[122].D45[sample])/(etatarray[214].D45[sample]+etatarray[214].R45[sample]-etatarray[122].D45[sample]-etatarray[122].R45[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR July-September 45-64: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D65[sample]-etatarray[122].D65[sample])/(etatarray[214].D65[sample]+etatarray[214].R65[sample]-etatarray[122].D65[sample]-etatarray[122].R65[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR July-September 65-74: mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[214].D75H[sample]-etatarray[122].D75H[sample])/(etatarray[214].D75H[sample]+etatarray[214].R75[sample]+etatarray[MAX_DAYS_DATA].RH[sample]-etatarray[122].D75H[sample]-etatarray[122].R75[sample]-etatarray[122].RH[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"IFR July-September 75+ (nurs. homes included): mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));



for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[24].R[sample])/(NPOP-etatarray[24].D[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune on March 30 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[39].R[sample])/(NPOP-etatarray[39].D[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune on April 14 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[MAX_DAYS_DATA].R[sample])/(NPOP-etatarray[MAX_DAYS_DATA].D[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune on day : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[24].AIR[sample])/(etatarray[24].S25[sample] + etatarray[24].S45[sample] + etatarray[24].S65[sample] + etatarray[24].E25[sample] + etatarray[24].E45[sample] + etatarray[24].E65[sample] + etatarray[24].AI25[sample] + etatarray[24].AI45[sample] + etatarray[24].AI65[sample]  + etatarray[24].PI25[sample] + etatarray[24].PI45[sample] + etatarray[24].PI65[sample] + etatarray[24].AIR[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune assympt on March 30 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));
printf("Perc immune assympt on March 30 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[39].AIR[sample])/(etatarray[39].S25[sample] + etatarray[39].S45[sample] + etatarray[39].S65[sample] + etatarray[39].E25[sample] + etatarray[39].E45[sample] + etatarray[39].E65[sample] + etatarray[39].AI25[sample] + etatarray[39].AI45[sample] + etatarray[39].AI65[sample]  + etatarray[39].PI25[sample] + etatarray[39].PI45[sample] + etatarray[39].PI65[sample] + etatarray[39].AIR[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune assympt on April 14 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));
printf("Perc immune assympt on April 14 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[MAX_DAYS_DATA-7].AIR[sample])/(etatarray[MAX_DAYS_DATA-7].S25[sample] + etatarray[MAX_DAYS_DATA-7].S45[sample] + etatarray[MAX_DAYS_DATA-7].S65[sample] + etatarray[MAX_DAYS_DATA-7].E25[sample] + etatarray[MAX_DAYS_DATA-7].E45[sample] + etatarray[MAX_DAYS_DATA-7].E65[sample] + etatarray[MAX_DAYS_DATA-7].AI25[sample] + etatarray[MAX_DAYS_DATA-7].AI45[sample] + etatarray[MAX_DAYS_DATA-7].AI65[sample]  + etatarray[MAX_DAYS_DATA-7].PI25[sample] + etatarray[MAX_DAYS_DATA-7].PI45[sample] + etatarray[MAX_DAYS_DATA-7].PI65[sample] + etatarray[MAX_DAYS_DATA-7].AIR[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune assympt on day : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95)); 

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[24].RH[sample])/(NPOPH);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune in homes on March 30 (vs total pop homes) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[39].RH[sample])/(NPOPH);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune in homes on April 14 (vs total pop homes) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[MAX_DAYS_DATA].RH[sample])/(NPOPH);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune in homes on day (vs total pop homes) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


	fprintf(result_par, "\n");


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[307].R[sample])/(NPOP-etatarray[307].D[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune on January 1 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95)); 


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[307].AIR[sample])/(etatarray[307].S25[sample] + etatarray[307].S45[sample] + etatarray[307].S65[sample] + etatarray[307].E25[sample] + etatarray[307].E45[sample] + etatarray[307].E65[sample] + etatarray[307].AI25[sample] + etatarray[307].AI45[sample] + etatarray[307].AI65[sample]  + etatarray[307].PI25[sample] + etatarray[307].PI45[sample] + etatarray[307].PI65[sample] + etatarray[307].AIR[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune assympt on January 1 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[307].RH[sample])/(NPOPH);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune in homes on January 1 (vs total pop homes) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[397].R[sample])/(NPOP-etatarray[397].D[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune on April 1 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95)); 



fprintf(result_par,"\n\n");

	for (looppar = 1; looppar <= NUMPAR; looppar++) {
	    tempsarray = foreachpararray(looppar, &pararray);
	    fprintf(result_par, "%s:      %.12f & %.12f & [%.12f ; %.12f]  \n",
		    foreachparname(looppar),
		    gsl_stats_mean(tempsarray, 1, SIZE_SAMPLE),
		    gsl_stats_median_from_sorted_data(tempsarray, 1, SIZE_SAMPLE),
		    gsl_stats_quantile_from_sorted_data(tempsarray, 1, SIZE_SAMPLE, 0.05), gsl_stats_quantile_from_sorted_data(tempsarray, 1, SIZE_SAMPLE, 0.95));
	}


	fclose(result_par);
	
	
		FILE *reimpout = NULL;
	reimpout = fopen("reimpout.txt", "w+");
	
	 for (day = DAY_START; day <= DAY_MAX; day++) {
	  

		fprintf(reimpout, "%d ", day);
		
for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= (NPOPY+NPOP25+NPOP45+NPOP65) * parsamples[sample].reimp * reimpfct(day,reimp)/NPOP;
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(reimpout," mean=%.1f med=%.1f [%.1f ; %.1f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

		}



    fclose(reimpout);
    
 


#if PER_CHANGE > 0 
    }
    }
#endif

return;

}