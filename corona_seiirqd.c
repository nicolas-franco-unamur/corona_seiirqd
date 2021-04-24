/*
Corona seiirqd covid-19 programm
Nicolas Franco - UNamur
preliminary undocumented version
*/


// day1 = 01/03   day32 = 01/04
#define DAY_START 1
#define DAY_MAX 397 //490		// 397- 490
#define MAX_DAYS_DATA 244     // 244
#define JUST_COMPUTE 0    // if 1 only perform MH with standard output (preval to just plot)
#define JUST_PLOT_MCMC 1   // if 1 prepare for plot from file
#define OUTPUT_FILE "out.txt"
#define SPEED   1 // 5
#define MOYRAND 1

// to simulate end of lockdown
#define PER_CHANGE_FROM 0
#define PER_CHANGE_TO 0.01
#define PER_CHANGE_STEP 10
#define PER_CHANGE2_FROM 1
#define PER_CHANGE2_TO  1.01
#define PER_CHANGE2_STEP 1
#define PER_CHANGE 1
#define WITH_ARG 0
#define PARTPLOT 1
#define PARTNAME "measures_quick"
// priors
#define READ_PRIORS 1
#define NUM_PRIORS 2900
// burn-in phase
#define NUM_BURNIN 20000	//  (100000 max 2h speed 5 - 10000)  h
#define SIZE_SAMPLE_BURNIN 2500	// 5 (1)
#define MAX_BURNIN 0
#define OPTIMIZE 5   // 5
// sample phase
#define NUM_ITE 20000		//  (2500-5000)
#define SIZE_SAMPLE_PER_BURNIN 1	// 10 (1)
#define BEST_FIT 0
// coef
#define COEF_SGM_BURNIN 5	//  5-10  (opt 100-200)
#define COEF_SGM_START_BURNIN 0	// 50 - 100
#define COEF_SGM_START 0	// 0  (0)
#define SIZE_SAMPLE (SIZE_SAMPLE_BURNIN*SIZE_SAMPLE_PER_BURNIN)
#define MALUS_LIKELIHOOD 10000000


#include "corona.h"

// compute likelihood

double likelihood(SEIIRQD * etat, SEIIRQD * etat_incid, SEIIRQD * data, SEIIRQD * data_incid, parameters * par)
{

    double likelihood = 0;
    int day;

    for (day = DAY_START; day <= MAX_DAYS_DATA; day++) {
	if (data[day].SI > 0 && data[day].SI > (NPOP-etat[day].S)) {
	    return 0.0;
	}

	
	// hospitals
	
	
	if (data_incid[day].SIQ > 0 && etat_incid[day].SIQ > 0) {  // with corrected data
	likelihood +=etat_incid[day].SIQ - (data_incid[day].SIQ*par->supphosp) * log(etat_incid[day].SIQ);
	 }

	
	if (data_incid[day].QR > 0 && etat_incid[day].QR > 0) {
	 likelihood +=etat_incid[day].QR - data_incid[day].QR * log(etat_incid[day].QR);
	}
	if (data_incid[day].QD > 0 && etat_incid[day].QD > 0) {
	 likelihood +=etat_incid[day].QD - data_incid[day].QD * log(etat_incid[day].QD);

	}
	
	
	// deaths


	if (data_incid[day].DY > 0  && etat_incid[day].DY > 0) {
	 likelihood +=etat_incid[day].DY - data_incid[day].DY * log(etat_incid[day].DY);
	}
	if (data_incid[day].D25 > 0  && etat_incid[day].D25 > 0) {
	 likelihood +=etat_incid[day].D25 - data_incid[day].D25 * log(etat_incid[day].D25);

	}	
	if (data_incid[day].D45 > 0  && etat_incid[day].D45 > 0) {
	 likelihood +=etat_incid[day].D45 - data_incid[day].D45 * log(etat_incid[day].D45);
	}
	if (data_incid[day].D65 > 0  && etat_incid[day].D65 > 0) {
	 likelihood +=etat_incid[day].D65 - data_incid[day].D65 * log(etat_incid[day].D65);

	}
	if (data_incid[day].D75 > 0  && etat_incid[day].D75 > 0) {
	 likelihood +=etat_incid[day].D75 - data_incid[day].D75 * log(etat_incid[day].D75);

	}
	
	if (data_incid[day].DH > 0  && etat_incid[day].DH > 0) {
	 likelihood +=etat_incid[day].DH - data_incid[day].DH * log(etat_incid[day].DH);

	}
	
	if (data_incid[day].D > 0  && etat_incid[day].D > 0) {
	 likelihood +=etat_incid[day].D - data_incid[day].D * log(etat_incid[day].D);
 
	}

	
    }

// tests

    // 30 march  (7 days delay)
    double pcinffromAI;
    pcinffromAI = (etat[24].AIR) / (etat[24].S25 + etat[24].S45 + etat[24].S65 + etat[24].E25 +etat[24].E45 + etat[24].E65 + etat[24].AI25 + etat[24].AI45 + etat[24].AI65 + etat[24].PI25 + etat[24].PI45 + etat[24].PI65 + etat[24].AIR);
    if (pcinffromAI < 0.005 || pcinffromAI > 0.028) {
	likelihood += MALUS_LIKELIHOOD * sqrt(pow(pcinffromAI - 0.013, 2));
    }

    // 14 april  (7 days delay)
    pcinffromAI = (etat[39].AIR) / (etat[39].S25 + etat[39].S45 + etat[39].S65 + etat[39].E25 +etat[39].E45 + etat[39].E65 + etat[39].AI25 + etat[39].AI45 + etat[39].AI65 + etat[39].PI25 + etat[39].PI45 + etat[39].PI65 + etat[39].AIR);
    if (pcinffromAI < 0.035|| pcinffromAI > 0.062) {  
	likelihood += MALUS_LIKELIHOOD * sqrt(pow(pcinffromAI - 0.047, 2));
    }
    
    // 27 april  (7 days delay)
    pcinffromAI = (etat[52].AIR) / (etat[52].S25 + etat[52].S45 + etat[52].S65 + etat[52].E25 +etat[52].E45 + etat[52].E65 + etat[52].AI25 + etat[52].AI45 + etat[52].AI65 + etat[52].PI25 + etat[52].PI45 + etat[52].PI65 + etat[52].AIR);
    if (pcinffromAI < 0.035|| pcinffromAI > 0.062) {  
	likelihood += MALUS_LIKELIHOOD * sqrt(pow(pcinffromAI - 0.047, 2));
    }

    // 11 mai  (7 days delay)
    pcinffromAI = (etat[66].AIR) / (etat[66].S25 + etat[66].S45 + etat[66].S65 + etat[66].E25 +etat[66].E45 + etat[66].E65 + etat[66].AI25 + etat[66].AI45 + etat[66].AI65 + etat[66].PI25 + etat[66].PI45 + etat[66].PI65 + etat[66].AIR);
    if (pcinffromAI < 0.035|| pcinffromAI > 0.062) {  
	likelihood += MALUS_LIKELIHOOD * sqrt(pow(pcinffromAI - 0.047, 2));
    }
    
     
    // rest homes asymptomatic
    double pcassh, ucount = 0, bcount = 0;
    for (day = 46; day <= MAX_DAYS_DATA; day++) {
	ucount += etat[day].EH + etat[day].AIH + etat[day].PIH;
	bcount += etat[day].EH + etat[day].AIH + etat[day].PIH + etat[day].SIH;
    }
    pcassh = ucount / bcount;
    if (pcassh < 0.75 - 0.10 || pcassh > 0.75 + 0.10) {
    likelihood += MALUS_LIKELIHOOD * sqrt(pow(pcassh - 0.75, 2));
    }


    double infectcount, popcount, infecth;
// rest homes infected 15-30 april  
    infectcount = 0;
    popcount = 0;
    for (day = 47; day < 63; day++) {
	infectcount += etat[day].EH + etat[day].AIH + etat[day].PIH + etat[day].SIH;
	popcount += etat[day].SH + etat[day].EH + etat[day].AIH + etat[day].PIH + etat[day].SIH + etat[day].RH;
    }
    infecth = infectcount / popcount;
    if (infecth < 0.08 - 0.03 || infecth > 0.08 + 0.03) {
    likelihood += MALUS_LIKELIHOOD * sqrt(pow(infecth - 0.08, 2));
    }


    // rest homes infected 15-31 May   
    infectcount = 0;
    popcount = 0;
    for (day = 77; day < 92; day++) {
	infectcount += etat[day].EH + etat[day].AIH + etat[day].PIH + etat[day].SIH;
	popcount += etat[day].SH + etat[day].EH + etat[day].AIH + etat[day].PIH + etat[day].SIH + etat[day].RH;
    }
    infecth = infectcount / popcount;
    if (infecth > 0.02 + 0.02) {
    likelihood += MALUS_LIKELIHOOD * sqrt(pow(infecth - 0.02, 2));
    }


    return likelihood;

}

// simu continue

int func(double t, const double y[], double f[], void *params)
{
    parameters *par = (parameters *) params;
    // y
    f[0] = -y[0]*par->lambdaa*(par->Myy*(y[2]+y[3])/NPOPY+par->My25*(y[9]+y[10])/NPOP25+par->My45*(y[16]+y[17])/NPOP45+par->My65*(y[23]+y[24])/NPOP65+par->My75*(y[30]+y[31])/NPOP75)
     -y[0]*par->lambdas*(par->Myy*y[4]/NPOPY+par->My25*y[11]/NPOP25+par->My45*y[18]/NPOP45+par->My65*y[25]/NPOP65+par->My75*y[32]/NPOP75);
    f[1] = y[0]*par->lambdaa*(par->Myy*(y[2]+y[3])/NPOPY+par->My25*(y[9]+y[10])/NPOP25+par->My45*(y[16]+y[17])/NPOP45+par->My65*(y[23]+y[24])/NPOP65+par->My75*(y[30]+y[31])/NPOP75)
     +y[0]*par->lambdas*(par->Myy*y[4]/NPOPY+par->My25*y[11]/NPOP25+par->My45*y[18]/NPOP45+par->My65*y[25]/NPOP65+par->My75*y[32]/NPOP75) - par->sigma * y[1];
    f[2] = par->pay *par->sigma * y[1] - par->gammaay * y[2];
    f[3] = (1-par->pay)*par->sigma * y[1] - par->tau * y[3];
    f[4] = par->tau * y[3] - par->deltay * y[4] - par->gammasy * y[4];
    f[5] = par->deltay * y[4] - par->gammaqyused * y[5] - par->ryused * y[5];
    f[6] = par->ryused * y[5];
    // 25
    f[7] = -y[7]*par->lambdaa*(par->M25y*(y[2]+y[3])/NPOPY+par->M2525*(y[9]+y[10])/NPOP25+par->M2545*(y[16]+y[17])/NPOP45+par->M2565*(y[23]+y[24])/NPOP65+par->M2575*(y[30]+y[31])/NPOP75)
     -y[7]*par->lambdas*(par->M25y*y[4]/NPOPY+par->M2525*y[11]/NPOP25+par->M2545*y[18]/NPOP45+par->M2565*y[25]/NPOP65+par->M2575*y[32]/NPOP75);
    f[8] = y[7]*par->lambdaa*(par->M25y*(y[2]+y[3])/NPOPY+par->M2525*(y[9]+y[10])/NPOP25+par->M2545*(y[16]+y[17])/NPOP45+par->M2565*(y[23]+y[24])/NPOP65+par->M2575*(y[30]+y[31])/NPOP75)
     +y[7]*par->lambdas*(par->M25y*y[4]/NPOPY+par->M2525*y[11]/NPOP25+par->M2545*y[18]/NPOP45+par->M2565*y[25]/NPOP65+par->M2575*y[32]/NPOP75) - par->sigma * y[8];
    f[9] = par->pa25 *par->sigma * y[8] - par->gammaa25 * y[9];
    f[10] = (1-par->pa25)*par->sigma * y[8] - par->tau * y[10];
    f[11] = par->tau * y[10] - par->delta25 * y[11] - par->gammas25 * y[11];
    f[12] = par->delta25 * y[11] - par->gammaq25used * y[12] - par->r25used * y[12];
    f[13] = par->r25used * y[12];  
    // 45
    f[14] = -y[14]*par->lambdaa*(par->M45y*(y[2]+y[3])/NPOPY+par->M4525*(y[9]+y[10])/NPOP25+par->M4545*(y[16]+y[17])/NPOP45+par->M4565*(y[23]+y[24])/NPOP65+par->M4575*(y[30]+y[31])/NPOP75)
     -y[14]*par->lambdas*(par->M45y*y[4]/NPOPY+par->M4525*y[11]/NPOP25+par->M4545*y[18]/NPOP45+par->M4565*y[25]/NPOP65+par->M4575*y[32]/NPOP75);
    f[15] = y[14]*par->lambdaa*(par->M45y*(y[2]+y[3])/NPOPY+par->M4525*(y[9]+y[10])/NPOP25+par->M4545*(y[16]+y[17])/NPOP45+par->M4565*(y[23]+y[24])/NPOP65+par->M4575*(y[30]+y[31])/NPOP75)
     +y[14]*par->lambdas*(par->M45y*y[4]/NPOPY+par->M4525*y[11]/NPOP25+par->M4545*y[18]/NPOP45+par->M4565*y[25]/NPOP65+par->M4575*y[32]/NPOP75) - par->sigma * y[15];
    f[16] = par->pa45 *par->sigma * y[15] - par->gammaa45 * y[16];
    f[17] = (1-par->pa45)*par->sigma * y[15] - par->tau * y[17];
    f[18] = par->tau * y[17] - par->delta45 * y[18] - par->gammas45 * y[18];
    f[19] = par->delta45 * y[18] - par->gammaq45used * y[19] - par->r45used * y[19];
    f[20] = par->r45used * y[19]; 
    // 65
    f[21] = -y[21]*par->lambdaa*(par->M65y*(y[2]+y[3])/NPOPY+par->M6525*(y[9]+y[10])/NPOP25+par->M6545*(y[16]+y[17])/NPOP45+par->M6565*(y[23]+y[24])/NPOP65+par->M6575*(y[30]+y[31])/NPOP75)
     -y[21]*par->lambdas*(par->M65y*y[4]/NPOPY+par->M6525*y[11]/NPOP25+par->M6545*y[18]/NPOP45+par->M6565*y[25]/NPOP65+par->M6575*y[32]/NPOP75);
    f[22] = y[21]*par->lambdaa*(par->M65y*(y[2]+y[3])/NPOPY+par->M6525*(y[9]+y[10])/NPOP25+par->M6545*(y[16]+y[17])/NPOP45+par->M6565*(y[23]+y[24])/NPOP65+par->M6575*(y[30]+y[31])/NPOP75)
     +y[21]*par->lambdas*(par->M65y*y[4]/NPOPY+par->M6525*y[11]/NPOP25+par->M6545*y[18]/NPOP45+par->M6565*y[25]/NPOP65+par->M6575*y[32]/NPOP75) - par->sigma * y[22];
    f[23] = par->pa65 *par->sigma * y[22] - par->gammaa65 * y[23];
    f[24] = (1-par->pa65)*par->sigma * y[22] - par->tau * y[24];
    f[25] = par->tau * y[24] - par->delta65 * y[25] - par->gammas65 * y[25];
    f[26] = par->delta65 * y[25] - par->gammaq65used * y[26] - par->r65used * y[26];
    f[27] = par->r65used * y[26];  
    // 75
    f[28] = -y[28]*par->lambdaa*(par->M75y*(y[2]+y[3])/NPOPY+par->M7525*(y[9]+y[10])/NPOP25+par->M7545*(y[16]+y[17])/NPOP45+par->M7565*(y[23]+y[24])/NPOP65+par->M7575*(y[30]+y[31])/NPOP75)
     -y[28]*par->lambdas*(par->M75y*y[4]/NPOPY+par->M7525*y[11]/NPOP25+par->M7545*y[18]/NPOP45+par->M7565*y[25]/NPOP65+par->M7575*y[32]/NPOP75);
    f[29] = y[28]*par->lambdaa*(par->M75y*(y[2]+y[3])/NPOPY+par->M7525*(y[9]+y[10])/NPOP25+par->M7545*(y[16]+y[17])/NPOP45+par->M7565*(y[23]+y[24])/NPOP65+par->M7575*(y[30]+y[31])/NPOP75)
     +y[28]*par->lambdas*(par->M75y*y[4]/NPOPY+par->M7525*y[11]/NPOP25+par->M7545*y[18]/NPOP45+par->M7565*y[25]/NPOP65+par->M7575*y[32]/NPOP75) - par->sigma * y[29];
    f[30] = par->pa75 *par->sigma * y[29] - par->gammaa75 * y[30];
    f[31] = (1-par->pa75)*par->sigma * y[29] - par->tau * y[31];
    f[32] = par->tau * y[31] - par->delta75 * y[32] - par->gammas75 * y[32];
    f[33] = par->delta75 * y[32] - par->gammaq75used * y[33] - par->r75used * y[33];
    f[34] = par->r75used * y[33];    
    // additional count
    f[35] = par->gammaa25 * y[9]+par->gammaa45 * y[16]+par->gammaa65 * y[23];
    f[36] = par->gammas25 * y[11]+par->gammas45 * y[18]+par->gammas65 * y[25];
    f[37] = par->gammaqyused * y[5]+par->gammaq25used * y[12]+par->gammaq45used * y[19]+par->gammaq65used * y[26]+par->gammaq75used * y[33];
    f[38] = par->deltay * y[4]+par->delta25 * y[11]+par->delta45 * y[18]+par->delta65 * y[25]+par->delta75 * y[32];
    f[39] = par->ryused * y[5]+par->r25used * y[12]+par->r45used * y[19]+par->r65used * y[26]+par->r75used * y[33];
    return GSL_SUCCESS;
}


int funchome(double t, const double y[], double f[], void *params)
{
    parameters *par = (parameters *) params;   
    f[0] = -y[0]*par->lambdaa*par->mh*(y[3]+y[4])/NPOPPERH -y[0]*par->lambdas*par->mh*y[5]/NPOPPERH - par->rhtused * (1 - par->pconf) * y[4];
    f[1] = y[0]*par->lambdaa*par->mh*(y[3]+y[4])/NPOPPERH +y[0]*par->lambdas*par->mh*y[5]/NPOPPERH - par->sigma * y[1];
    f[2] = par->pah *par->sigma * y[1] - par->gammaah * y[2];
    f[3] = (1-par->pah)*par->sigma * y[1] - par->tau * y[3];
    f[4] = par->tau * y[3] - par->deltahused * y[4] - par->gammash * y[4] - par->rhtused * par->pconf * y[4];
    f[5] = par->deltahused * y[4] - par->gammaqhused * y[5] - par->rhused * y[5];
    f[6] = par->rhused * y[5];
    f[7] = par->rhtused * y[4];
    f[8] = par->deltahused * y[4];
    f[9] = par->gammaqhused * y[5];  
    f[10] = par->gammaah * y[2] + par->gammash * y[4] + par->gammaqhused * y[5];  
    return GSL_SUCCESS;
}

void contact_matrices(parameters * par, double work, double school, double leisure, double home){
    // home
    par->Myy = 1.5277237*home;
    par->My25 = 1.2068022*home;
    par->My45 = 0.9249841*home;
    par->My65 = 0.1132896*home;
    par->My75 = 0.08835537*home;
    par->M25y = 1.2922885*home;
    par->M2525 = 0.9838020*home;
    par->M2545 = 0.6479268*home;
    par->M2565 = 0.1702118*home;
    par->M2575 = 0.22757490*home;
    par->M45y = 0.9893866*home;
    par->M4525 = 0.6471937*home;
    par->M4545 = 1.1973751*home;
    par->M4565 = 0.1468694*home;
    par->M4575 = 0.23303846*home;
    par->M65y = 0.3920099*home;
    par->M6525 = 0.5500133*home;
    par->M6545 = 0.4751236*home;
    par->M6565 = 0.5834762*home;
    par->M6575 = 0.26491845*home;
    par->M75y = 0.2892186*home;
    par->M7525 = 0.6956556*home;
    par->M7545 = 0.7131637*home;
    par->M7565 = 0.2506101*home;
    par->M7575 = 1.00546883*home;
    // work (+transport)    
    par->Myy += work*0.47389136;
    par->My25 += work*0.5960882;
    par->My45 +=  work*0.3694173;
    par->My65 += work*0.02262112;
    par->My75 += work*0.006680393;
    par->M25y += work*0.63831329;
    par->M2525 += work*2.3903069;
    par->M2545 += work*1.6179651;
    par->M2565 += work*0.11526048;
    par->M2575 += work*0.090571246;
    par->M45y += work*0.39513818;
    par->M4525 += work*1.6161345;
    par->M4545 += work*1.4517584;
    par->M4565 += work*0.11615251; 
    par->M4575 += work*0.105063074;
    par->M65y += work*0.07827463;
    par->M6525 += work*0.3724466;
    par->M6545 += work*0.3757542;
    par->M6565 += work*0.27208083;
    par->M6575 += work*0.165808622;
    par->M75y += work*0.02186731;
    par->M7525 += work*0.2768600;
    par->M7545 += work*0.3215227;
    par->M7565 += work*0.15685324;
    par->M7575 += work*0.100006623;
    // school    
    par->Myy += school*5.03765528;
    par->My25 += school*0.34386299;
    par->My45 +=  school*0.16926068;
    par->My65 += school*0.004074279;
    par->My75 += school*0.000000000;
    par->M25y += school*0.36822123;
    par->M2525 += school*0.09900565;
    par->M2545 += school*0.03918760;
    par->M2565 += school*0.000000000; 
    par->M2575 += school*0.000000000;
    par->M45y += school*0.18104553;
    par->M4525 += school*0.03914326;
    par->M4545 += school*0.02996797;
    par->M4565 += school*0.006859308;
    par->M4575 += school*0.000000000;
    par->M65y += school*0.01409801;
    par->M6525 += school*0.00000000;
    par->M6545 += school*0.02218991;
    par->M6565 += school*0.000000000; 
    par->M6575 += school*0.000000000;
    par->M75y += school*0.00000000;
    par->M7525 += school*0.00000000;
    par->M7545 += school*0.00000000;
    par->M7565 += school*0.00000000;
    par->M7575 += school*0.000000000;
    // leisure (+otherplaces)   
    par->Myy += leisure*3.0841839;
    par->My25 += leisure*1.1718123;
    par->My45 +=  leisure*0.6688232;
    par->My65 += leisure*0.2125313;
    par->My75 += leisure*0.1264132;
    par->M25y += leisure*1.2548200;
    par->M2525 += leisure*2.8294124;
    par->M2545 += leisure*1.7006531;
    par->M2565 += leisure*0.3163111;
    par->M2575 += leisure*0.1747568;
    par->M45y += leisure*0.7153903;
    par->M4525 += leisure*1.6987289;
    par->M4545 += leisure*2.1215761;
    par->M4565 += leisure*0.4566617;
    par->M4575 += leisure*0.3233819;
    par->M65y += leisure*0.7354104;
    par->M6525 += leisure*1.0221109;
    par->M6545 += leisure*1.4773036;
    par->M6565 += leisure*1.1055254;
    par->M6575 += leisure*0.5253612;
    par->M75y += leisure*0.4137955;
    par->M7525 += leisure*0.5342003;
    par->M7545 += leisure*0.9896401;
    par->M7565 += leisure*0.4969863;
    par->M7575 += leisure*0.3103922;
}

double compute_rzero(parameters * par, double work, double school, double leisure, double home){
    contact_matrices(par, work, school, leisure,home);

 



double dataai[] = { par->pay*par->Myy/par->gammaay, par->pa25*par->My25/par->gammaa25, par->pa45*par->My45/par->gammaa45, par->pa65*par->My65/par->gammaa65, par->pa75*par->My75/par->gammaa75,
par->pay*par->M25y/par->gammaay, par->pa25*par->M2525/par->gammaa25, par->pa45*par->M2545/par->gammaa45, par->pa65*par->M2565/par->gammaa65, par->pa75*par->M2575/par->gammaa75,
par->pay*par->M45y/par->gammaay, par->pa25*par->M4525/par->gammaa25, par->pa45*par->M4545/par->gammaa45, par->pa65*par->M4565/par->gammaa65, par->pa75*par->M4575/par->gammaa75,
par->pay*par->M65y/par->gammaay, par->pa25*par->M6525/par->gammaa25, par->pa45*par->M6545/par->gammaa45, par->pa65*par->M6565/par->gammaa65, par->pa75*par->M6575/par->gammaa75,
par->pay*par->M75y/par->gammaay, par->pa25*par->M7525/par->gammaa25, par->pa45*par->M7545/par->gammaa45, par->pa65*par->M7565/par->gammaa65, par->pa75*par->M7575/par->gammaa75   };

  gsl_matrix_view mai   = gsl_matrix_view_array (dataai, 5, 5);
    
    gsl_matrix_scale(&mai.matrix, par->lambdaa);
    
double datapi[] = { (1-par->pay)*par->Myy/par->tau, (1-par->pa25)*par->My25/par->tau, (1-par->pa45)*par->My45/par->tau, (1-par->pa65)*par->My65/par->tau, (1-par->pa75)*par->My75/par->tau,
(1-par->pay)*par->M25y/par->tau, (1-par->pa25)*par->M2525/par->tau, (1-par->pa45)*par->M2545/par->tau, (1-par->pa65)*par->M2565/par->tau, (1-par->pa75)*par->M2575/par->tau,
(1-par->pay)*par->M45y/par->tau, (1-par->pa25)*par->M4525/par->tau, (1-par->pa45)*par->M4545/par->tau, (1-par->pa65)*par->M4565/par->tau, (1-par->pa75)*par->M4575/par->tau,
(1-par->pay)*par->M65y/par->tau, (1-par->pa25)*par->M6525/par->tau, (1-par->pa45)*par->M6545/par->tau, (1-par->pa65)*par->M6565/par->tau, (1-par->pa75)*par->M6575/par->tau,
(1-par->pay)*par->M75y/par->tau, (1-par->pa25)*par->M7525/par->tau, (1-par->pa45)*par->M7545/par->tau, (1-par->pa65)*par->M7565/par->tau, (1-par->pa75)*par->M7575/par->tau   };

  gsl_matrix_view mpi   = gsl_matrix_view_array (datapi, 5, 5);
   gsl_matrix_scale(&mpi.matrix, par->lambdaa);


double datasi[] = { (1-par->pay)*par->Myy/(par->gammasy+ par->deltay), (1-par->pa25)*par->My25/(par->gammas25+ par->delta25), (1-par->pa45)*par->My45/(par->gammas45+ par->delta45), (1-par->pa65)*par->My65/(par->gammas65+ par->delta65), (1-par->pa75)*par->My75/(par->gammas75+ par->delta75),
(1-par->pay)*par->M25y/(par->gammasy+ par->deltay), (1-par->pa25)*par->M2525/(par->gammas25+ par->delta25), (1-par->pa45)*par->M2545/(par->gammas45+ par->delta45), (1-par->pa65)*par->M2565/(par->gammas65+ par->delta65), (1-par->pa75)*par->M2575/(par->gammas75+ par->delta75),
(1-par->pay)*par->M45y/(par->gammasy+ par->deltay), (1-par->pa25)*par->M4525/(par->gammas25+ par->delta25), (1-par->pa45)*par->M4545/(par->gammas45+ par->delta45), (1-par->pa65)*par->M4565/(par->gammas65+ par->delta65), (1-par->pa75)*par->M4575/(par->gammas75+ par->delta75),
(1-par->pay)*par->M65y/(par->gammasy+ par->deltay), (1-par->pa25)*par->M6525/(par->gammas25+ par->delta25), (1-par->pa45)*par->M6545/(par->gammas45+ par->delta45), (1-par->pa65)*par->M6565/(par->gammas65+ par->delta65), (1-par->pa75)*par->M6575/(par->gammas75+ par->delta75),
(1-par->pay)*par->M75y/(par->gammasy+ par->deltay), (1-par->pa25)*par->M7525/(par->gammas25+ par->delta25), (1-par->pa45)*par->M7545/(par->gammas45+ par->delta45), (1-par->pa65)*par->M7565/(par->gammas65+ par->delta65), (1-par->pa75)*par->M7575/(par->gammas75+ par->delta75)   };


  gsl_matrix_view msi   = gsl_matrix_view_array (datasi, 5, 5);
    
    gsl_matrix_scale(&msi.matrix, par->lambdas);
    
   gsl_matrix_add(&msi.matrix, &mpi.matrix);
   gsl_matrix_add(&msi.matrix, &mai.matrix);

  gsl_vector_complex *eval = gsl_vector_complex_alloc (5);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (5, 5);

  gsl_eigen_nonsymmv_workspace * w =
    gsl_eigen_nonsymmv_alloc (5);

  gsl_eigen_nonsymmv (&msi.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_DESC);
  
  double return_value;
                           
  return_value= GSL_REAL(gsl_vector_complex_get (eval, 0));
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
  
  return return_value;
                           
}

double reimpfct(int day, double reimp[4][DAY_MAX + 1]){
    int delay=0;
    
    if(day<124){
        return 0;
    }       
    // July
    if(day<155){
        return 0.36*(0.23*reimp[1][day-delay]+0.11*reimp[2][day-delay]+0.09*reimp[3][day-delay]+0.07*reimp[4][day-delay]);
    } 
    // August
    if(day<186){
        return 0.26*(0.23*reimp[1][day-delay]+0.11*reimp[2][day-delay]+0.09*reimp[3][day-delay]+0.07*reimp[4][day-delay]);
    }  
    // September
    if(day<201){
        return 2.0/3.0*0.21*(0.23*reimp[1][day-delay]+0.11*reimp[2][day-delay]+0.09*reimp[3][day-delay]+0.07*reimp[4][day-delay]);
    }         
    if(day<216){
        return 1.0/3.0*0.21*(0.23*reimp[1][day-delay]+0.11*reimp[2][day-delay]+0.09*reimp[3][day-delay]+0.07*reimp[4][day-delay]);
    }          
    
    
    return 0;
}

void perform_simu(int start, int end, int openhomes, double *removedfromgenpop,SEIIRQD * etat, SEIIRQD * etat_incid, parameters * par, gsl_rng * r, double *y, double z[NHOMES_SIMU][NHOMESVAR], gsl_odeiv2_driver * d, gsl_odeiv2_driver * dhome, double reimp[4][DAY_MAX + 1])
{
    int day = 1, nhome,loop;
    double t,*temp_incid,*temp_daypo,*temp_day, recovery_coefp=1,recovery_coefm=1;


    for (day = start; day < end; day++) {
	t = day;
	
	
	y[1] += NPOPY * par->reimp * reimpfct(day,reimp)/NPOP;
    y[8] += NPOP25 * par->reimp * reimpfct(day,reimp)/NPOP;
    y[15] += NPOP45 * par->reimp * reimpfct(day,reimp)/NPOP;
    y[22] += NPOP65 * par->reimp * reimpfct(day,reimp)/NPOP;
    y[0] -= NPOPY * par->reimp * reimpfct(day,reimp)/NPOP;
    y[7] -= NPOP25 * par->reimp * reimpfct(day,reimp)/NPOP;
    y[14] -= NPOP45 * par->reimp * reimpfct(day,reimp)/NPOP;
    y[21] -= NPOP65 * par->reimp * reimpfct(day,reimp)/NPOP;
		
	recovery_coefp = 1 + par->recovery / (1 + exp(-(day - par->rcap) / par->rcaps));
	recovery_coefm = 1 - par->recovery / (1 + exp(-(day - par->rcap) / par->rcaps));
	

	par->gammaqyused=par->gammaqy*recovery_coefp;
	par->gammaq25used=par->gammaq25*recovery_coefp;
	par->gammaq45used=par->gammaq45*recovery_coefp;
	par->gammaq65used=par->gammaq65*recovery_coefp;
	par->gammaq75used=par->gammaq75*recovery_coefp;
	par->gammaqhused=par->gammaqh*recovery_coefp;
    par->ryused=par->ry*recovery_coefm;
    par->r25used=par->r25*recovery_coefm;
    par->r45used=par->r45*recovery_coefm;
    par->r65used=par->r65*recovery_coefm;
    par->r75used=par->r75*recovery_coefm;
    par->rhused=par->rh*recovery_coefm;
	
	if(day-par->delay<1){
	par->deltahused = par->deltah;
	par->rhtused = 0; 
	} else if(day<124){
	par->deltahused = par->deltah - par->rht * par->pconf / (1 + exp(-(etat[day-(int)par->delay].Q - par->cap) / par->caps));
	par->rhtused = par->rht / (1 + exp(-(etat[day-(int)par->delay].Q - par->cap) / par->caps));
	} else {
	   	par->deltahused = par->deltah - par->rht * par->pconf / (1 + exp(-(0 - par->cap) / par->caps));
	par->rhtused = par->rht / (1 + exp(-(0 - par->cap) / par->caps));

	}
	
	gsl_odeiv2_driver_apply(d, &t, day + 1, y);
	for (nhome = 0; nhome < NHOMES_SIMU; nhome++) {
	    t = day;
	    gsl_odeiv2_driver_apply(dhome, &t, day + 1, z[nhome]);
	}
	etat[day + 1].SY = y[0];
	etat[day + 1].EY = y[1];
	etat[day + 1].AIY = y[2];
	etat[day + 1].PIY = y[3];
	etat[day + 1].SIY = y[4];
	etat[day + 1].QY = y[5];
	etat[day + 1].DY = y[6];
    etat[day + 1].S25 = y[7];
	etat[day + 1].E25 = y[8];
	etat[day + 1].AI25 = y[9];
	etat[day + 1].PI25 = y[10];
	etat[day + 1].SI25 = y[11];
	etat[day + 1].Q25 = y[12];
	etat[day + 1].D25 = y[13];
	etat[day + 1].S45 = y[14];
	etat[day + 1].E45 = y[15];
	etat[day + 1].AI45 = y[16];
	etat[day + 1].PI45 = y[17];
	etat[day + 1].SI45 = y[18];
	etat[day + 1].Q45 = y[19];
	etat[day + 1].D45 = y[20];
	etat[day + 1].S65 = y[21];
	etat[day + 1].E65 = y[22];
	etat[day + 1].AI65 = y[23];
	etat[day + 1].PI65 = y[24];
	etat[day + 1].SI65 = y[25];
	etat[day + 1].Q65 = y[26];
	etat[day + 1].D65 = y[27];
	etat[day + 1].S75 = y[28];
	etat[day + 1].E75 = y[29];
	etat[day + 1].AI75 = y[30];
	etat[day + 1].PI75 = y[31];
	etat[day + 1].SI75 = y[32];
	etat[day + 1].Q75 = y[33];
	etat[day + 1].D75 = y[34];
	etat[day + 1].AIR = y[35];
	etat[day + 1].SIR = y[36];
    etat[day + 1].QR = y[37];
	etat[day + 1].SIQ = y[38];
	etat[day + 1].QD = y[39];
	etat[day + 1].SH = 0;
	etat[day + 1].EH = 0;
	etat[day + 1].AIH = 0;
	etat[day + 1].SIH = 0;
	etat[day + 1].PIH = 0;
	etat[day + 1].QH = 0;
	etat[day + 1].DH = 0;
	etat[day + 1].RH = 0;
	etat[day + 1].HQ = 0;
	for (nhome = 0; nhome < NHOMES_SIMU; nhome++) {
	    z[nhome][0] = z[nhome][0] > 0 ? z[nhome][0] : 0;
	    z[nhome][1] = z[nhome][1] > 0 ? z[nhome][1] : 0;
	    z[nhome][2] = z[nhome][2] > 0 ? z[nhome][2] : 0;
	    z[nhome][3] = z[nhome][3] > 0 ? z[nhome][3] : 0;
	    z[nhome][4] = z[nhome][4] > 0 ? z[nhome][4] : 0;
	    z[nhome][5] = z[nhome][5] > 0 ? z[nhome][5] : 0;
// add infected with proba
        if(openhomes==1){
	         if ((par->pth / NPOP ) * z[nhome][0] *  
            (par->lambdaa* (etat[day].AIY + etat[day].PIY +etat[day].AI25 + etat[day].PI25 + etat[day].AI45 + etat[day].PI45 + etat[day].AI65 + etat[day].PI65 + etat[day].AI75 + etat[day].PI75)
            +par->lambdas* (etat[day].SIY + etat[day].SI25 + etat[day].SI45 + etat[day].SI65 + etat[day].SI75) )
             > gsl_rng_uniform(r)) {  
		       if (z[nhome][0] >= 1) {
		        z[nhome][0]--;
		        z[nhome][1]++;
		       }
	        }
        } else {
            if ((par->pthp / NPOP ) * z[nhome][0] *  
            (par->lambdaa* (etat[day].AI25 + etat[day].PI25 + etat[day].AI45 + etat[day].PI45)
            +par->lambdas* (etat[day].SI25 + etat[day].SI45) )
             > gsl_rng_uniform(r)) {
		      if (z[nhome][0] >= 1) {
		        z[nhome][0]--;
		        z[nhome][1]++;
		      }
	       }
        }
	    if (y[28] > 10000) {
		*removedfromgenpop += ((NPOPPERH - z[nhome][1] - z[nhome][2] - z[nhome][3] - z[nhome][4] - z[nhome][5] - z[nhome][10]) - z[nhome][0])*COEFSPEED;
		y[28] -= ((NPOPPERH - z[nhome][1] - z[nhome][2] - z[nhome][3] - z[nhome][4] - z[nhome][5] - z[nhome][10]) - z[nhome][0])*COEFSPEED;
		z[nhome][0] = NPOPPERH - z[nhome][1] - z[nhome][2] - z[nhome][3] - z[nhome][4] - z[nhome][5] - z[nhome][10];
	    }
	    etat[day + 1].SH += z[nhome][0]*COEFSPEED;
	    etat[day + 1].EH += z[nhome][1]*COEFSPEED;
	    etat[day + 1].AIH += z[nhome][2]*COEFSPEED;
	    etat[day + 1].PIH += z[nhome][3]*COEFSPEED;
	    etat[day + 1].SIH += z[nhome][4]*COEFSPEED;
	    etat[day + 1].QH += z[nhome][5]*COEFSPEED;
	    etat[day + 1].D75 += z[nhome][6]*COEFSPEED;
	    etat[day + 1].QD += z[nhome][6]*COEFSPEED;
	    etat[day + 1].DH += z[nhome][7]*COEFSPEED;
	    etat[day + 1].SIQ += z[nhome][8]*COEFSPEED;
	    etat[day + 1].HQ += z[nhome][8]*COEFSPEED;
	    etat[day + 1].QR += z[nhome][9]*COEFSPEED;
	    etat[day + 1].RH += z[nhome][10]*COEFSPEED;
	   	}
	etat[day + 1].RY = NPOPY - etat[day + 1].SY - etat[day + 1].EY - etat[day + 1].AIY - etat[day + 1].SIY - etat[day + 1].PIY - etat[day + 1].QY - etat[day + 1].DY;
	etat[day + 1].R25 = NPOP25 - etat[day + 1].S25 - etat[day + 1].E25 - etat[day + 1].AI25 - etat[day + 1].SI25 - etat[day + 1].PI25 - etat[day + 1].Q25 - etat[day + 1].D25;
	etat[day + 1].R45 = NPOP45 - etat[day + 1].S45 - etat[day + 1].E45 - etat[day + 1].AI45 - etat[day + 1].SI45 - etat[day + 1].PI45 - etat[day + 1].Q45 - etat[day + 1].D45;
	etat[day + 1].R65 = NPOP65 - etat[day + 1].S65 - etat[day + 1].E65 - etat[day + 1].AI65 - etat[day + 1].SI65 - etat[day + 1].PI65 - etat[day + 1].Q65 - etat[day + 1].D65;
	etat[day + 1].S75 = y[28];
	etat[day + 1].R75 = NPOP75 - *removedfromgenpop - etat[day + 1].S75 - etat[day + 1].E75 - etat[day + 1].AI75 - etat[day + 1].SI75  - etat[day + 1].PI75 - etat[day + 1].Q75 - y[34];
	etat[day + 1].S = etat[day + 1].SY + etat[day + 1].S25 + etat[day + 1].S45 + etat[day + 1].S65 + etat[day + 1].S75 + etat[day + 1].SH;
	etat[day + 1].E = etat[day + 1].EY + etat[day + 1].E25 + etat[day + 1].E45 + etat[day + 1].E65 + etat[day + 1].E75 + etat[day + 1].EH;
	etat[day + 1].AI = etat[day + 1].AIY + etat[day + 1].AI25 + etat[day + 1].AI45 + etat[day + 1].AI65 + etat[day + 1].AI75 + etat[day + 1].AIH;
	etat[day + 1].PI = etat[day + 1].PIY + etat[day + 1].PI25 + etat[day + 1].PI45 + etat[day + 1].PI65 + etat[day + 1].PI75 + etat[day + 1].PIH;
	etat[day + 1].SI = etat[day + 1].SIY + etat[day + 1].SI25 + etat[day + 1].SI45 + etat[day + 1].SI65 + etat[day + 1].SI75 + etat[day + 1].SIH;
	etat[day + 1].R = etat[day + 1].RY + etat[day + 1].R25 + etat[day + 1].R45 + etat[day + 1].R65 + etat[day + 1].R75 + etat[day + 1].RH;
	etat[day + 1].Q = etat[day + 1].QY + etat[day + 1].Q25 + etat[day + 1].Q45 + etat[day + 1].Q65 + etat[day + 1].Q75 + etat[day + 1].QH;
	etat[day + 1].D = etat[day + 1].DY + etat[day + 1].D25 + etat[day + 1].D45 + etat[day + 1].D65 + etat[day + 1].D75 + etat[day + 1].DH;
	etat[day + 1].D75H = etat[day + 1].D75 + etat[day + 1].DH;
	
// incidence
    for (loop = 1; loop <= NUMSEG; loop++) {
	temp_day = foreachetat(loop, &etat[day]);
	temp_daypo = foreachetat(loop, &etat[day+1]);
	temp_incid = foreachetat(loop, &etat_incid[day+1]);
	*temp_incid = *temp_daypo-*temp_day;
    }
    
       
    }
    

}


void simucontinu(SEIIRQD * etat, SEIIRQD * etat_incid, parameters * par, int end, gsl_rng * r, double reimp[4][DAY_MAX + 1])
{

    int  nhome;
    double  removedfromgenpop = 0;
    double y[40];
    double z[NHOMES_SIMU][NHOMESVAR];

    y[0] = etat[DAY_START].SY;
    y[1] = etat[DAY_START].EY;
    y[2] = etat[DAY_START].AIY;
    y[3] = etat[DAY_START].PIY;
    y[4] = etat[DAY_START].SIY;
    y[5] = etat[DAY_START].QY;
    y[6] = etat[DAY_START].DY;
    y[7] = etat[DAY_START].S25;
    y[8] = etat[DAY_START].E25;
    y[9] = etat[DAY_START].AI25;
    y[10] = etat[DAY_START].PI25;
    y[11] = etat[DAY_START].SI25;
    y[12] = etat[DAY_START].Q25;
    y[13] = etat[DAY_START].D25;    
    y[14] = etat[DAY_START].S45;
    y[15] = etat[DAY_START].E45;
    y[16] = etat[DAY_START].AI45;
    y[17] = etat[DAY_START].PI45;
    y[18] = etat[DAY_START].SI45;
    y[19] = etat[DAY_START].Q45;
    y[20] = etat[DAY_START].D45;    
    y[21] = etat[DAY_START].S65;
    y[22] = etat[DAY_START].E65;
    y[23] = etat[DAY_START].AI65;
    y[24] = etat[DAY_START].PI65;
    y[25] = etat[DAY_START].SI65;
    y[26] = etat[DAY_START].Q65;
    y[27] = etat[DAY_START].D65;   
    y[28] = etat[DAY_START].S75;
    y[29] = etat[DAY_START].E75;
    y[30] = etat[DAY_START].AI75;
    y[31] = etat[DAY_START].PI75;
    y[32] = etat[DAY_START].SI75;
    y[33] = etat[DAY_START].Q75;
    y[34] = etat[DAY_START].D75;  
    y[35] = etat[DAY_START].AIR;
    y[36] = etat[DAY_START].SIR;
    y[37] = etat[DAY_START].QR;
    y[38] = etat[DAY_START].SIQ;
    y[39] = etat[DAY_START].QD;

    for (nhome = 0; nhome < NHOMES_SIMU; nhome++) {
	z[nhome][0] = NPOPPERH;
	z[nhome][1] = 0;
	z[nhome][2] = 0;
	z[nhome][3] = 0;
	z[nhome][4] = 0;
	z[nhome][5] = 0;
	z[nhome][6] = 0;
	z[nhome][7] = 0;
	z[nhome][8] = 0;
	z[nhome][9] = 0;
	z[nhome][10] = 0;
    }

    contact_matrices(par, 1, 1, 1,1);

    gsl_odeiv2_system sys = { func, NULL, 40, par };
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-1, 1e-1, 0.0);
    gsl_odeiv2_system syshome = { funchome, NULL, 11, par };
    gsl_odeiv2_driver *dhome = gsl_odeiv2_driver_alloc_y_new(&syshome, gsl_odeiv2_step_rkf45, 1e-0, 1e-1, 0.0);

    perform_simu(DAY_START,15,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);

// lockdown
     contact_matrices(par, 1, 0, par->leisure_lock,1);
     perform_simu(15,20,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_lock, 0, par->leisure_lock,par->home_lock);
    perform_simu(20,66,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
 
// lift of lock
    contact_matrices(par, (par->work_lock+par->work_unlock)/2, 0, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(66,73,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, 0, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(73,80,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, 0.2*par->school_unlock, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(80,87,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, 0.4*par->school_unlock, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(87,95,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp); 
    contact_matrices(par, par->work_unlock, 0.6*par->school_unlock, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(95,101,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, par->school_unlock, par->leisure_june,par->home_unlock);
    perform_simu(101,124,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // July
    contact_matrices(par, 0.2*par->work_unlock, 0,par->leisure_july,par->home_lock+(par->home_unlock-par->home_lock)*2);  // modif work and home
    perform_simu(124,152,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // August
    contact_matrices(par, 0.2*par->work_unlock, 0,par->leisure_augustus,par->home_unlock);
    perform_simu(152,186,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);      
    // september
    contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(186,210,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_lock+(par->home_unlock-par->home_lock)*2); // buble relaxed
    perform_simu(210,224,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
// first new measures october
        contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(224,248,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // toussain
        contact_matrices(par, par->work_unlock, 0,par->leisure_sept,par->home_unlock);
    perform_simu(248,255,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
      //
        contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(255,297,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // christmas
        contact_matrices(par, par->work_unlock, 0,par->leisure_sept,par->home_unlock);
    perform_simu(297,311,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // january
        contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(311,353,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // carnaval
            contact_matrices(par, par->work_unlock, 0,par->leisure_sept,par->home_unlock);
    perform_simu(353,360,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    //
        contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(360,399,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);

            contact_matrices(par, par->work_unlock, 0,par->leisure_sept,par->home_unlock);
    perform_simu(399,415,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
        contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(415,end,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    


       
    gsl_odeiv2_driver_free(d);
    gsl_odeiv2_driver_free(dhome);

}

// for simu end of lockdown
void simucontinu2(SEIIRQD * etat, SEIIRQD * etat_incid, parameters * par, int end, double change,double change2, gsl_rng * r, double reimp[4][DAY_MAX + 1])
{


   int  nhome;
    double  removedfromgenpop = 0;
    double y[40];
    double z[NHOMES_SIMU][NHOMESVAR];

    y[0] = etat[DAY_START].SY;
    y[1] = etat[DAY_START].EY;
    y[2] = etat[DAY_START].AIY;
    y[3] = etat[DAY_START].PIY;
    y[4] = etat[DAY_START].SIY;
    y[5] = etat[DAY_START].QY;
    y[6] = etat[DAY_START].DY;
    y[7] = etat[DAY_START].S25;
    y[8] = etat[DAY_START].E25;
    y[9] = etat[DAY_START].AI25;
    y[10] = etat[DAY_START].PI25;
    y[11] = etat[DAY_START].SI25;
    y[12] = etat[DAY_START].Q25;
    y[13] = etat[DAY_START].D25;    
    y[14] = etat[DAY_START].S45;
    y[15] = etat[DAY_START].E45;
    y[16] = etat[DAY_START].AI45;
    y[17] = etat[DAY_START].PI45;
    y[18] = etat[DAY_START].SI45;
    y[19] = etat[DAY_START].Q45;
    y[20] = etat[DAY_START].D45;    
    y[21] = etat[DAY_START].S65;
    y[22] = etat[DAY_START].E65;
    y[23] = etat[DAY_START].AI65;
    y[24] = etat[DAY_START].PI65;
    y[25] = etat[DAY_START].SI65;
    y[26] = etat[DAY_START].Q65;
    y[27] = etat[DAY_START].D65;   
    y[28] = etat[DAY_START].S75;
    y[29] = etat[DAY_START].E75;
    y[30] = etat[DAY_START].AI75;
    y[31] = etat[DAY_START].PI75;
    y[32] = etat[DAY_START].SI75;
    y[33] = etat[DAY_START].Q75;
    y[34] = etat[DAY_START].D75;  
    y[35] = etat[DAY_START].AIR;
    y[36] = etat[DAY_START].SIR;
    y[37] = etat[DAY_START].QR;
    y[38] = etat[DAY_START].SIQ;
    y[39] = etat[DAY_START].QD;

    for (nhome = 0; nhome < NHOMES_SIMU; nhome++) {
	z[nhome][0] = NPOPPERH;
	z[nhome][1] = 0;
	z[nhome][2] = 0;
	z[nhome][3] = 0;
	z[nhome][4] = 0;
	z[nhome][5] = 0;
	z[nhome][6] = 0;
	z[nhome][7] = 0;
	z[nhome][8] = 0;
	z[nhome][9] = 0;
	z[nhome][10] = 0;
    }

    contact_matrices(par, 1, 1, 1,1);

    gsl_odeiv2_system sys = { func, NULL, 40, par };
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-1, 1e-1, 0.0);
    gsl_odeiv2_system syshome = { funchome, NULL, 11, par };
    gsl_odeiv2_driver *dhome = gsl_odeiv2_driver_alloc_y_new(&syshome, gsl_odeiv2_step_rkf45, 1e-0, 1e-1, 0.0);

    perform_simu(DAY_START,15,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);

// lockdown
     contact_matrices(par, 1, 0, par->leisure_lock,1);
     perform_simu(15,20,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_lock, 0, par->leisure_lock,par->home_lock);
    perform_simu(20,66,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
 
// lift of lock
      contact_matrices(par, (par->work_lock+par->work_unlock)/2, 0, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(66,73,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, 0, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(73,80,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, 0.2*par->school_unlock, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(80,87,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, 0.4*par->school_unlock, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(87,95,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp); 
    contact_matrices(par, par->work_unlock, 0.6*par->school_unlock, par->leisure_lock,par->home_lock+(par->home_unlock-par->home_lock)/2);
    perform_simu(95,101,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, par->school_unlock, par->leisure_june,par->home_unlock);
    perform_simu(101,124,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // July
    contact_matrices(par, 0.2*par->work_unlock, 0,par->leisure_july,par->home_lock+(par->home_unlock-par->home_lock)*2);  // modif work and home
    perform_simu(124,152,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // August
    contact_matrices(par, 0.2*par->work_unlock, 0,par->leisure_augustus,par->home_unlock);
    perform_simu(152,186,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);      
    // september
    contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(186,210,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_lock+(par->home_unlock-par->home_lock)*2); // buble relaxed
    perform_simu(210,224,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
       
 // double newschool =  par->school_sept; 
  double newwork = par->work_unlock;
  double newleisure = par->leisure_sept;
  double newhome = par->home_unlock;
  //newschool = 5.0/6.0*par->school_sept;
  newwork = par->work_lock + (par->work_unlock - par->work_lock)*2/4;
  newleisure = par->leisure_lock+(par->leisure_sept-par->leisure_lock)*2/4;
  newhome = par->home_lock + (par->home_unlock - par->home_lock)*2/4;
  //newwork = par->work_lock;
  //newleisure = par->leisure_lock;
  //newhome = par->home_lock;

  double longschool =  par->school_sept; 
  double longwork = par->work_unlock;
  double longleisure = par->leisure_sept;
  double longhome = par->home_unlock;
  longwork = par->work_lock + (par->work_unlock - par->work_lock)*0/4;
  longleisure = par->leisure_lock+(par->leisure_sept-par->leisure_lock)*0/4;
  longhome = par->home_lock + (par->home_unlock - par->home_lock)*0/4;
  longschool = 5.0/6.0*par->school_sept;
 // longwork = par->work_lock;
  //longleisure = par->leisure_lock;
  //longhome = par->home_lock;
  
    // (first new measures october)
        contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(224,234,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
// new measures 19 october
      contact_matrices(par, newwork, par->school_sept,newleisure,newhome);
    perform_simu(234,248,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);    
        // toussain (extended)
        contact_matrices(par, longwork, 0,longleisure,longhome);
    perform_simu(248,262,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
     // school reopened
             contact_matrices(par, longwork, longschool,longleisure,longhome);
    perform_simu(262,290,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
   // end of measures
    contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(290,297,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // christmas 
    contact_matrices(par, par->work_unlock, 0,par->leisure_sept,par->home_unlock);
    perform_simu(297,311,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    
      double scschool =  par->school_sept; 
  double scwork = par->work_unlock;
  double scleisure = par->leisure_sept;
  double schome = par->home_unlock;

//scleisure +=change/100.0;
 // scschool=1;  
 //schome +=change/100.0;
  // scwork +=change/100.0;
    
    // january (restart)
    contact_matrices(par, scwork, scschool,scleisure,schome);
    perform_simu(311,353,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    // carnaval
    contact_matrices(par, scwork, 0,scleisure,schome);
    perform_simu(353,360,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    //
    contact_matrices(par, scwork, scschool,scleisure,schome);
    perform_simu(360,399,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);

           contact_matrices(par, scwork, 0,scleisure,schome);
    perform_simu(399,415,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
        contact_matrices(par, scwork, scschool,scleisure,schome);
    perform_simu(415,end,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    

       
    gsl_odeiv2_driver_free(d);
    gsl_odeiv2_driver_free(dhome);  
}



void init(SEIIRQD * etat, SEIIRQD * etat_incid, parameters * par)
{
    int loop;
    double *temp;
    for (loop = 1; loop <= NUMSEG; loop++) {
	temp = foreachetat(loop, &etat[DAY_START]);
	*temp = 0;
	temp = foreachetat(loop, &etat_incid[DAY_START]);
	*temp = 0;
    }
    etat[DAY_START].EY = NPOPY * par->p0;
    etat[DAY_START].E25 = NPOP25 * par->p0;
    etat[DAY_START].E45 = NPOP45 * par->p0;
    etat[DAY_START].E65 = NPOP65 * par->p0;
    etat[DAY_START].E75 = NPOP75 * par->p0;
    etat[DAY_START].E = etat[DAY_START].EY + etat[DAY_START].E25 + etat[DAY_START].E45 + etat[DAY_START].E65 + etat[DAY_START].E75;
    etat[DAY_START].SY = NPOPY - etat[DAY_START].EY;
    etat[DAY_START].S25 = NPOP25 - etat[DAY_START].E25;
    etat[DAY_START].S45 = NPOP45 - etat[DAY_START].E45;
    etat[DAY_START].S65 = NPOP65 - etat[DAY_START].E65;
    etat[DAY_START].S75 = NPOP75 - etat[DAY_START].E75;
    etat[DAY_START].SH = NPOPH;
    etat[DAY_START].S = NPOP - etat[DAY_START].E;

}

double simulikelihood(SEIIRQD * etat, SEIIRQD * etat_incid, SEIIRQD * data, SEIIRQD * data_incid, parameters * par, gsl_rng * r, double reimp[4][DAY_MAX + 1]){
   double templikelihood=0;

    
    if(MOYRAND>1){
        int i;
        templikelihood=0;
        for(i=0;i<MOYRAND;i++){
            init(etat, etat_incid, par);
    	    simucontinu(etat, etat_incid, par, MAX_DAYS_DATA, r,reimp);
            templikelihood +=likelihood(etat, etat_incid, data, data_incid,par);	
        }   
        templikelihood = templikelihood/MOYRAND;
        #if	JUST_COMPUTE == 0
		    fprintf(stderr, " %.0f\r", templikelihood);
#endif
        return templikelihood;
    }
    

    	init(etat, etat_incid, par);
    	simucontinu(etat, etat_incid, par, MAX_DAYS_DATA, r,reimp);
	templikelihood=likelihood(etat, etat_incid, data, data_incid,par);	
	
        #if	JUST_COMPUTE == 0
		    fprintf(stderr, " %.0f  \r", templikelihood);
#endif	
	return templikelihood;
		
}


void readdata(SEIIRQD * data, SEIIRQD * data_incid, double reimp[4][DAY_MAX + 1])
{
    FILE *datafile = NULL,*reimpfile = NULL;
    int day, readday;
    double totalD75,totalD75_incid;



    datafile = fopen("data.txt", "r");
    if (datafile == NULL) {
	printf("Fail read data\n ");
	exit(-1);
    }
    for (day = 1; day <= MAX_DAYS_DATA; day++) {
	fscanf(datafile, "%d", &readday);
	if (!(readday == day)) {
	    printf("Error read day %d-%d\n ", readday, day);
	    exit(-1);
	}
	fscanf(datafile, "%lf", &data[day].SI);
	fscanf(datafile, "%lf %lf",&data_incid[day].SIQ, &data[day].SIQ);
    fscanf(datafile, "%lf %lf",&data_incid[day].Q, &data[day].Q);
	fscanf(datafile, "%lf %lf",&data_incid[day].QR, &data[day].QR);
	fscanf(datafile, "%lf %lf",&data_incid[day].QD, &data[day].QD);
	fscanf(datafile, "%lf %lf",&data_incid[day].DH, &data[day].DH);
	fscanf(datafile, "%lf %lf",&data_incid[day].D, &data[day].D);
	fscanf(datafile, "%lf %lf %lf %lf %lf",&data_incid[day].DY, &data_incid[day].D25, &data_incid[day].D45, &data_incid[day].D65, &totalD75_incid);
	fscanf(datafile, "%lf %lf %lf %lf %lf",&data[day].DY, &data[day].D25, &data[day].D45, &data[day].D65, &totalD75);
	data_incid[day].D75 = totalD75_incid - data_incid[day].DH;
	data[day].D75 = totalD75 - data[day].DH;
    }

    fclose(datafile);
    
    for (day = 1; day <= DAY_MAX; day++) {
       reimp[1][day] = 0;
       reimp[2][day] = 0;
       reimp[3][day] = 0;
       reimp[4][day] = 0;
    }    
    

    
    reimpfile = fopen("reimp.txt", "r");
    if (reimpfile == NULL) {
	printf("Fail read data\n ");
	exit(-1);
    }
    for(day=MAX_DAYS_DATA; day>=124; day-- ){
    fscanf(reimpfile, "%d", &day);
    fscanf(reimpfile, "%lf", &reimp[1][day]);
    fscanf(reimpfile, "%lf", &reimp[2][day]);
    fscanf(reimpfile, "%lf", &reimp[3][day]);
    fscanf(reimpfile, "%lf", &reimp[4][day]);
   
    } 
    
    fclose(reimpfile);
    

}



// fct test par
int test_par(parameters par){
if(par.p0*NPOP <1 || par.p0*NPOP>20000){return 0;}
if(par.lambdaa <=0){return 0;}
if(par.lambdas <par.lambdaa){return 0;}
if(par.sigma <0.25 || par.sigma>1){return 0;}
if(par.tau <0.1 || par.tau>0.5){return 0;}
if(par.pay <par.pa25 || par.pay>1 ){return 0;}
if(par.pa25 <par.pa45){return 0;}
if(par.pa45 <par.pa65){return 0;}
if(par.pa65 <par.pa75){return 0;}
if(par.pa75 <par.pah){return 0;}
if(par.pah <0){return 0;}
if(par.deltay <0 || par.deltay>0.5){return 0;}
if(par.delta25 < par.deltay || par.delta25>0.6){return 0;}
if(par.delta45 < par.delta25 || par.delta45>0.6){return 0;}
if(par.delta65 < par.delta45 || par.delta65>0.7){return 0;}
if(par.delta75 < par.delta65 || par.delta75>0.8){return 0;}
if(par.deltah < par.delta75 || par.deltah>0.8){return 0;}
if(par.gammaay <0.1 || par.gammaay>0.5){return 0;}
if(par.gammaa25 <0 || par.gammaa25>par.gammaay){return 0;}
if(par.gammaa45 <0 || par.gammaa45>par.gammaa25){return 0;}
if(par.gammaa65 <0 || par.gammaa65>par.gammaa45){return 0;}
if(par.gammaa75 <0 || par.gammaa75>par.gammaa65){return 0;}
if(par.gammaah <0 || par.gammaah>par.gammaa75){return 0;}
if(par.gammasy <0.1 || par.gammasy>0.5){return 0;}
if(par.gammas25 <0 || par.gammas25>par.gammasy){return 0;}
if(par.gammas45 <0 || par.gammas45>par.gammas25){return 0;}
if(par.gammas65 <0 || par.gammas65>par.gammas45){return 0;}
if(par.gammas75 <0 || par.gammas75>par.gammas65){return 0;}
if(par.gammash <0 || par.gammash>par.gammas75){return 0;}
if(par.gammaqy <0 || par.gammaqy>par.gammasy){return 0;}
if(par.gammaq25 <0 || par.gammaq25>par.gammaqy){return 0;}
if(par.gammaq45 <0 || par.gammaq45>par.gammaq25){return 0;}
if(par.gammaq65 <0 || par.gammaq65>par.gammaq45){return 0;}
if(par.gammaq75 <0 || par.gammaq75>par.gammaq65){return 0;}
if(par.gammaqh <0 || par.gammaqh>par.gammaq75){return 0;}
if(par.ry <0 || par.ry>0.05){return 0;}
if(par.r25 <par.ry || par.r25>0.1){return 0;}
if(par.r45 <par.r25 || par.r45>0.1){return 0;}
if(par.r65 <par.r45 || par.r65>0.15){return 0;}
if(par.r75 <par.r65 || par.r75>0.2){return 0;}
if(par.rh <par.r75 || par.rh>0.5){return 0;}
if(par.rht <0 || par.rht*par.pconf>par.deltah){return 0;}
if(par.cap <0){return 0;}
if(par.caps <=0){return 0;}
if(par.pconf<0.5 || par.pconf>0.9){return 0;}
if(par.pthp<0){return 0;}
if(par.pth<par.pthp){return 0;}
if(par.mh<0){return 0;}
if(par.supphosp<1 || par.supphosp>1.3){return 0;}
if(par.delay<0 || par.delay>60){return 0;} 
if(par.home_lock <0.25 || par.home_lock >1){return 0;}
if(par.work_lock <0 || par.work_lock>1){return 0;}
if(par.leisure_lock <0 || par.leisure_lock>1){return 0;}
if(par.home_unlock <par.home_lock || par.home_unlock>1){return 0;}
if(par.work_unlock <par.work_lock || par.work_unlock>1){return 0;}
if(par.school_unlock <0 || par.school_unlock>0.5){return 0;}
if(par.leisure_june <par.leisure_lock || par.leisure_june>1){return 0;}
if(par.leisure_july <0 || par.leisure_july>1){return 0;} 
if(par.leisure_augustus<0 || par.leisure_augustus >par.leisure_july){return 0;} 
if(par.leisure_sept<par.leisure_augustus || par.leisure_sept >1){return 0;}
if(par.school_sept <par.school_unlock || par.school_sept>1){return 0;}
if(par.recovery <0 || par.recovery>0.9){return 0;}
if(par.rcap <0){return 0;}
if(par.rcaps <=0){return 0;}
if(par.reimp <0){return 0;}


double rzero,rzeropp;
rzero=compute_rzero(&par,1,1,1,1);
if(rzero <2.5 || rzero> 4.5){return 0;}
rzeropp=compute_rzero(&par,par.work_lock,0,par.leisure_lock,par.home_lock);
if(rzeropp <0.5 || rzeropp> 1){return 0;}
return 1;

}



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


parameters random_step_one(gsl_rng * r){
    double coef=1;
    parameters par;
    int looppar,num;
    
    for (looppar = 1; looppar <= NUMPAR; looppar++) {
		*foreachpar(looppar, &par)=0;		
	}
	
	num=gsl_rng_get(r) % NUMPAR;
	


    
    
	
/*	int force=gsl_rng_get(r)%20;  // attention temps
	
	switch(force){
	  case 0: 	num=60; break;  // attention temps 
	  case 1: 	num=61; break;  // attention temps 
	  case 2: 	num=62; break;  // attention temps 
	  case 3: 	num=63; break;  // attention temps 
	  case 4: 	num=64; break;  // attention temps 
	}*/
		

        
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

parameters opti_randomize(parameters parburninsamples, SEIIRQD * data, SEIIRQD * data_incid, double *burninlikelihood, int burnin, gsl_rng * r,double reimp[4][DAY_MAX + 1]){
    parameters  parbn,par_step,testpar;
	SEIIRQD etatbn[DAY_MAX + 1], etatbn_incid[DAY_MAX + 1];
	int accept,looppar,i;
	double newlikelihoodbn=0,testlikelihoodbn=0,coef;
	// new opti test reverse
	
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
		    fprintf(stderr, "              lklh... (opti) %.0f  (%.0f%%)  pth%f  pthp%f  R0H%f  ju%f  j%f  a%f  s%f  reim%f  sch%f  \r", *burninlikelihood,100.0 * (burnin) / NUM_BURNIN,
		    parbn.pth,parbn.pthp,parbn.mh/(parbn.pah*parbn.gammaah +(1-parbn.pah)*(parbn.gammash+parbn.deltah)),
parbn.leisure_june,parbn.leisure_july,parbn.leisure_augustus,parbn.leisure_sept,parbn.reimp,parbn.school_sept);
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
		    fprintf(stderr, "              lklh..r (opti) %.0f  (%.0f%%)  pth%f  pthp%f  R0H%f  ju%f  j%f  a%f  s%f  reim%f  sch%f   \r",  *burninlikelihood,100.0 * (burnin) / NUM_BURNIN,
		    parbn.pth,parbn.pthp,parbn.mh/(parbn.pah*parbn.gammaah +(1-parbn.pah)*(parbn.gammash+parbn.deltah)),
parbn.leisure_june,parbn.leisure_july,parbn.leisure_augustus,parbn.leisure_sept,parbn.reimp,parbn.school_sept);
#endif


		    return parbn;		    
		} 
		}
    
    
        return parburninsamples;
}



unsigned long int random_seed()
{

    unsigned int seed;
    struct timeval tv;
    FILE *devrandom;

    if ((devrandom = fopen("/dev/random", "r")) == NULL) {
	gettimeofday(&tv, 0);
	seed = tv.tv_sec + tv.tv_usec;
    } else {
	fread(&seed, sizeof(seed), 1, devrandom);
	fclose(devrandom);
    }

    return (seed);

}

// *************************************************************************************

// fct main
int main(int argc, char **argv)
{
    SEIIRQD data[DAY_MAX + 1], data_incid[DAY_MAX + 1];
    double reimp[4][DAY_MAX + 1];
    int looppar;
    const gsl_rng_type *T;
    gsl_rng *r;
// init
    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, random_seed());
    readdata(data, data_incid,reimp);

//  Monte-Carlo Metropolis-Hastings *************************************************************************************
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
	printf("erreur taille specifiee priors: %d < %d\n", num_prior, NUM_PRIORS);
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


// start burning 
    for (samplebn = 0; samplebn < SIZE_SAMPLE_BURNIN; samplebn++) {

	int sample = 0, badburnin = 1, burnin = 0, accept = 0;
	double  burninlikelihood = 0, burninlikelihoodcopy = 0, bnlklh0 = 0;
	parameters parburninsamples;

	while (badburnin) {	// to be sure to be in an admissible zone
	   burninlikelihood = 0;
	   burninlikelihoodcopy = 0;
	   bnlklh0 = 0;

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


#if READ_PRIORS ==1
	    parburninsamples = priors_array[gsl_rng_get(r) % NUM_PRIORS]; 
#endif

burninlikelihood = DBL_MAX;

for (burnin = 0; burnin < NUM_BURNIN; burnin++) {

	
#if OPTIMIZE > 0	
    parburninsamples=opti_randomize(parburninsamples, data, data_incid, &burninlikelihood,burnin,r,reimp);
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

// end of burnin    ****************************************************
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

#if BEST_FIT == 1
		if (newlikelihood < finallikelihood) {	// for FC-HC random walk
#else

		if (newlikelihood < finallikelihood ||  gsl_rng_uniform(r) < exp(- newlikelihood + finallikelihood)) {	// for MH gsl_rng_uniform(r)
						   if(newlikelihood > finallikelihood ){ fprintf(stderr, "                lklh accept %.10f                \n", exp(- newlikelihood + finallikelihood));}  // temp


#endif

		    finallikelihood = newlikelihood;
		    temp_par = par;
		    bnlklh0 = likelihood(etat, etat_incid, data, data_incid,&par);
		    #if JUST_COMPUTE == 0
		    fprintf(stderr, "                lklh sample %.0f    (%.0f%%)               \r",bnlklh0,100.0 * (ite) / NUM_ITE);

#endif
		}

	    }			// fin for ite
// print results
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

#endif				//   end just plot 0 = end of MCMC-MH     *******************************************************************
#if JUST_COMPUTE ==0

// read parameters
    parameters parsamples[SIZE_SAMPLE];
    FILE *parameters_mcmc;
    int sample, num_sample, day;

    parameters_mcmc = fopen("inputparameters.txt", "r");
    fscanf(parameters_mcmc, "%d", &num_sample);
    if (num_sample < SIZE_SAMPLE) {
	printf("erreur taille specifiee: %d < %d\n", num_sample, SIZE_SAMPLE);
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

	   //	rzeroh[sample] =parsamples[sample].mh/(parsamples[sample].pah*parsamples[sample].gammaah+(1-parsamples[sample].pah)*(parsamples[sample].gammash+parsamples[sample].deltah));

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


// plot stats     ****************************************************
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




/*
for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tauy/(parsamples[sample].gammay+parsamples[sample].tauy)*parsamples[sample].deltay/(parsamples[sample].gammay+parsamples[sample].deltay);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of hospitalized 0-44 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tau45/(parsamples[sample].gamma45+parsamples[sample].tau45)*parsamples[sample].delta45/(parsamples[sample].gamma45+parsamples[sample].delta45);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of hospitalized 45-64 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tau65/(parsamples[sample].gamma65+parsamples[sample].tau65)*parsamples[sample].delta65/(parsamples[sample].gamma65+parsamples[sample].delta65);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of hospitalized 65-74 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tau75/(parsamples[sample].gamma75+parsamples[sample].tau75)*parsamples[sample].delta75/(parsamples[sample].gamma75+parsamples[sample].delta75);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of hospitalized 75+ : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tauh/(parsamples[sample].gammah+parsamples[sample].tauh)*parsamples[sample].deltah/(parsamples[sample].gammah+parsamples[sample].deltah);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of hospitalized Homes (max) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tauh/(parsamples[sample].gammah+parsamples[sample].tauh)*(parsamples[sample].deltah-parsamples[sample].rht*parsamples[sample].pconf)/(parsamples[sample].gammah+parsamples[sample].deltah);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of hospitalized Homes (min) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tauy/(parsamples[sample].gammay+parsamples[sample].tauy)*parsamples[sample].deltay/(parsamples[sample].gammay+parsamples[sample].deltay)*(parsamples[sample].ry)/(parsamples[sample].ry+parsamples[sample].epsilony);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality 0-44 : mean=%.4f%% med=%.4f%% [%.4f%% ; %.4f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tau45/(parsamples[sample].gamma45+parsamples[sample].tau45)*parsamples[sample].delta45/(parsamples[sample].gamma45+parsamples[sample].delta45)*(parsamples[sample].r45)/(parsamples[sample].r45+parsamples[sample].epsilon45);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality 45-64 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tau65/(parsamples[sample].gamma65+parsamples[sample].tau65)*parsamples[sample].delta65/(parsamples[sample].gamma65+parsamples[sample].delta65)*(parsamples[sample].r65)/(parsamples[sample].r65+parsamples[sample].epsilon65);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality 65-74 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tau75/(parsamples[sample].gamma75+parsamples[sample].tau75)*parsamples[sample].delta75/(parsamples[sample].gamma75+parsamples[sample].delta75)*(parsamples[sample].r75)/(parsamples[sample].r75+parsamples[sample].epsilon75);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality 75+ : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tauh/(parsamples[sample].gammah+parsamples[sample].tauh)*
((parsamples[sample].deltah)/(parsamples[sample].gammah+parsamples[sample].deltah)*(parsamples[sample].rh)/(parsamples[sample].rh+parsamples[sample].epsilonh));
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality Homes (March) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*parsamples[sample].tauh/(parsamples[sample].gammah+parsamples[sample].tauh)*
((parsamples[sample].deltah-parsamples[sample].rht*parsamples[sample].pconf)/(parsamples[sample].gammah+parsamples[sample].deltah)*(parsamples[sample].rh)/(parsamples[sample].rh+parsamples[sample].epsilonh) +  (parsamples[sample].rht*parsamples[sample].pconf/(parsamples[sample].gammah+parsamples[sample].deltah)));
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality Homes (April) : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));



for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(parsamples[sample].ry)/(parsamples[sample].ry+parsamples[sample].epsilony);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality in hospital 0-44 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(parsamples[sample].r45)/(parsamples[sample].r45+parsamples[sample].epsilon45);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality in hospital 45-64 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(parsamples[sample].r65)/(parsamples[sample].r65+parsamples[sample].epsilon65);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality in hospital 65-74 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(parsamples[sample].r75)/(parsamples[sample].r75+parsamples[sample].epsilon75);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality in hospital 75+ : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(parsamples[sample].rh)/(parsamples[sample].rh+parsamples[sample].epsilonh);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc of mortality in hospital homes : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));*/



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

/*for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= etatarray[DAY_MAX].D[sample];
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Deaths on end : mean=%.0f med=%.0f [%.0f ; %.0f]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*etatarray[DAY_MAX].D[sample]/(NPOP);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Deaths on end : mean=%.3f%% med=%.3f%% [%.3f%% ; %.3f%%]  \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));*/

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

/*#if PER_CHANGE == 1


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]= 100*(etatarray[DAY_MAX].R[sample])/(NPOP-etatarray[DAY_MAX].D[sample]);
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"Perc immune on 2021 : mean=%.2f%% med=%.2f%% [%.2f%% ; %.2f%%]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"\n");

fprintf(result_par,"\\item Effective $R_t$ June 8: mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotjune , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotjune , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotjune, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotjune , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"\\item Effective $R_t$ July 1: mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotjuly , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotjuly , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotjuly, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotjuly , 1, SIZE_SAMPLE,0.95));

fprintf(result_par,"\\item Effective $R_t$ Sept. 1: mean=%.6f med=%.6f [%.6f ; %.6f]  \n",
gsl_stats_mean(rzerotseptember , 1, SIZE_SAMPLE),
gsl_stats_median_from_sorted_data(rzerotseptember , 1, SIZE_SAMPLE),
gsl_stats_quantile_from_sorted_data(rzerotseptember, 1, SIZE_SAMPLE,0.05),
gsl_stats_quantile_from_sorted_data(rzerotseptember , 1, SIZE_SAMPLE,0.95));

double stat2[SIZE_SAMPLE];

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]=0;
stat2[sample]=0;
for(day=60; day<= DAY_MAX-15; day++){
if( etatarray[day-15].Q[sample]<etatarray[day].Q[sample] && etatarray[day-1].Q[sample]<etatarray[day].Q[sample]  && etatarray[day-2].Q[sample]<etatarray[day].Q[sample] && etatarray[day-3].Q[sample]<etatarray[day].Q[sample]  &&  etatarray[day].Q[sample]>etatarray[day+1].Q[sample] && etatarray[day].Q[sample]>etatarray[day+2].Q[sample] && etatarray[day].Q[sample]>etatarray[day+3].Q[sample] && etatarray[day].Q[sample]>etatarray[day+15].Q[sample] ){
stat[sample]=etatarray[day].Q[sample];
break;
}
}
for(day+=15; day<= DAY_MAX-15; day++){
if( etatarray[day-15].Q[sample]<etatarray[day].Q[sample] && etatarray[day-1].Q[sample]<etatarray[day].Q[sample]  && etatarray[day-2].Q[sample]<etatarray[day].Q[sample] && etatarray[day-3].Q[sample]<etatarray[day].Q[sample]  &&  etatarray[day].Q[sample]>etatarray[day+1].Q[sample] && etatarray[day].Q[sample]>etatarray[day+2].Q[sample] && etatarray[day].Q[sample]>etatarray[day+3].Q[sample] && etatarray[day].Q[sample]>etatarray[day+15].Q[sample] ){
stat2[sample]=etatarray[day].Q[sample];
break;
}
}

}

gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"\\item Second peak height: mean=%.0f med=%.0f [%.0f ; %.0f]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));


gsl_sort(stat2 ,1, SIZE_SAMPLE);
fprintf(result_par,"\\item Third peak height: mean=%.0f med=%.0f [%.0f ; %.0f]   \n",gsl_stats_mean(stat2 , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat2 , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat2, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat2 , 1, SIZE_SAMPLE,0.95));


for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]=etatarray[DAY_MAX].QR[sample]+etatarray[DAY_MAX].QD[sample] -(etatarray[93].QR[sample]+etatarray[93].QD[sample] );
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"\\item Total hospitalised (June1-Dec31):  mean={\\bf %.0f} med=%.0f [%.0f ; %.0f]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]=etatarray[DAY_MAX].QD[sample] -(etatarray[93].QD[sample] );
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"\\item Total deaths (June1-Dec31): mean=%.0f med=%.0f [%.0f ; %.0f]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]=etatarray[DAY_MAX].QR[sample] -(etatarray[93].QR[sample] );
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"\\item Total released (June1-Dec31): mean=%.0f med=%.0f [%.0f ; %.0f]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]=0;
for(day=94 ; day<= DAY_MAX ; day++){
stat[sample]+=etatarray[DAY_MAX].Q[sample];
}
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(result_par,"\\item Sum of numbers of patients each day (area of the curve, June1-Dec31):\\\\ mean=%.0f med=%.0f [%.0f ; %.0f]   \n",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));



#endif*/


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
    
    // output fagg
    
 /*  FILE * csv_output = NULL; 
    
    sprintf(namefile, "../report/fagg5/Unamur_Model_%s.csv",PARTNAME);
    csv_output =fopen(namefile,"w+");
    
fprintf(csv_output,"Date, incidences_mean, incidences_median, incidences_LL, incidences_UL, load_mean, load_median, load_LL, load_UL\n");


for(day=DAY_START+1;day<= DAY_MAX;day++){

fprintf(csv_output,"%d", day-1);

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]=etatarray_incid[day].SIQ[sample]/parsamples[sample].supphosp;
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(csv_output,",%.0f,%.0f,%.0f,%.0f",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.025),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.975));


tempsetatarray=etatarray[day].Q;
fprintf(csv_output,",%.0f,%.0f,%.0f,%.0f", gsl_stats_mean(tempsetatarray, 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(tempsetatarray, 1, SIZE_SAMPLE,0.025),gsl_stats_quantile_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE,0.975));

fprintf(csv_output,"\n");


}


    
    
fclose(csv_output);*/
    

// output Masquelier
/*#if PER_CHANGE > 0
FILE * csv_output = NULL;

perc=5;
percpouc = 0.01*(double)perc;
csv_output = NULL;

#if RELOCK ==1
sprintf (namefile, "deaths/deaths_R0%1.1f-1-%1.1f.csv",change, change );
#else
sprintf (namefile, "deaths/deaths_R0%1.1f.csv",change );
#endif

csv_output =fopen(namefile,"w+");

fprintf(csv_output,"\"Day (day1=March1)\"");
fprintf(csv_output,",\"Total deaths (mean)\",\"Total deaths (median)\",\"Total deaths (P5)\",\"Total deaths (P95)\"");
fprintf(csv_output,",\"0-44 (mean)\",\"0-44 (median)\",\"0-44 (P5)\",\"0-44 (P95)\"");
fprintf(csv_output,",\"45-64 (mean)\",\"45-64 (median)\",\"45-64 (P5)\",\"45-64 (P95)\"");
fprintf(csv_output,",\"65-74 (mean)\",\"65-74 (median)\",\"65-74 (P5)\",\"65-74 (P95)\"");
fprintf(csv_output,",\"75+ (mean)\",\"75+ (median)\",\"75+ (P5)\",\"75+ (P95)\"");
fprintf(csv_output,",\"Total deaths covid only (mean)\",\"Total deaths covid only (median)\",\"Total deaths covid only (P5)\",\"Total deaths covid only (P95)\"");
fprintf(csv_output,",\"75+ covid only (mean)\",\"75+ covid only (median)\",\"75+ covid only (P5)\",\"75+ covid only (P95)\"");
fprintf(csv_output,"\n");


for(day=DAY_START+1;day<= DAY_MAX;day++){

fprintf(csv_output,"%d", day-1);

tempsetatarray=etatarray[day].D;
fprintf(csv_output,",%.0f,%.0f,%.0f,%.0f", gsl_stats_mean(tempsetatarray, 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(tempsetatarray, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE,0.95));

tempsetatarray=etatarray[day].DY;
fprintf(csv_output,",%.0f,%.0f,%.0f,%.0f", gsl_stats_mean(tempsetatarray, 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(tempsetatarray, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE,0.95));

tempsetatarray=etatarray[day].D45;
fprintf(csv_output,",%.0f,%.0f,%.0f,%.0f", gsl_stats_mean(tempsetatarray, 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(tempsetatarray, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE,0.95));

tempsetatarray=etatarray[day].D65;
fprintf(csv_output,",%.0f,%.0f,%.0f,%.0f", gsl_stats_mean(tempsetatarray, 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(tempsetatarray, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE,0.95));

tempsetatarray=etatarray[day].D75H;
fprintf(csv_output,",%.0f,%.0f,%.0f,%.0f", gsl_stats_mean(tempsetatarray, 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(tempsetatarray, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(tempsetatarray , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]=etatarray[day].DH[sample]*parsamples[sample].pconf+etatarray[day].D75[sample]+etatarray[day].D65[sample]+etatarray[day].D45[sample]+etatarray[day].DY[sample];
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(csv_output,",%.0f,%.0f,%.0f,%.0f",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

for(sample=0; sample <SIZE_SAMPLE;sample++){
stat[sample]=etatarray[day].DH[sample]*parsamples[sample].pconf+etatarray[day].D75[sample];
}
gsl_sort(stat ,1, SIZE_SAMPLE);
fprintf(csv_output,",%.0f,%.0f,%.0f,%.0f",gsl_stats_mean(stat , 1, SIZE_SAMPLE),gsl_stats_median_from_sorted_data(stat , 1, SIZE_SAMPLE),gsl_stats_quantile_from_sorted_data(stat, 1, SIZE_SAMPLE,0.05),gsl_stats_quantile_from_sorted_data(stat , 1, SIZE_SAMPLE,0.95));

fprintf(csv_output,"\n");


}
fclose(csv_output);

#endif*/




#if PER_CHANGE > 0 
    }
    }
#endif


#endif				//   end just compute 0       *******************************************************************


    return 0;

}
