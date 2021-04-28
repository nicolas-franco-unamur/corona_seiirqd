/*
Corona seiirqd covid-19 program
File: simulation.c (function for simulation)
Copyright Nicolas Franco - UNamur 2020-2021
Public version 1.0 (28/04/2021)
Version corresponding to the scientific paper:
"Covid-19 Belgium: Extended SEIR-QD model with nursing homes and long-term scenarios-based forecasts"
https://doi.org/10.1101/2020.09.07.20190108
*/

/*
Technical functions for running the simulation
*/

#include "corona.h"

// compute likelihood

double likelihood(SEIIRQD * etat, SEIIRQD * etat_incid, SEIIRQD * data, SEIIRQD * data_incid, parameters * par)
{

    double likelihood = 0;
    int day;

    // tests part
    for (day = DAY_START; day <= MAX_DAYS_DATA; day++) {
	if (data[day].SI > 0 && data[day].SI > (NPOP-etat[day].S)) {
	    return 0.0;
	}

	// hospitals part
		
	if (data_incid[day].SIQ > 0 && etat_incid[day].SIQ > 0) {  // with corrected data
	likelihood +=etat_incid[day].SIQ - (data_incid[day].SIQ*par->supphosp) * log(etat_incid[day].SIQ);
	 }	
	if (data_incid[day].QR > 0 && etat_incid[day].QR > 0) {
	 likelihood +=etat_incid[day].QR - data_incid[day].QR * log(etat_incid[day].QR);
	}
	if (data_incid[day].QD > 0 && etat_incid[day].QD > 0) {
	 likelihood +=etat_incid[day].QD - data_incid[day].QD * log(etat_incid[day].QD);

	}

	// deaths part

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

// serological data

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
    
// nursing homes part     
    
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


// function for differential equations general population
int func(double t, const double y[], double f[], void *params)
{
    parameters *par = (parameters *) params;
    // y (0-24)
    f[0] = -y[0]*par->lambdaa*(par->Myy*(y[2]+y[3])/NPOPY+par->My25*(y[9]+y[10])/NPOP25+par->My45*(y[16]+y[17])/NPOP45+par->My65*(y[23]+y[24])/NPOP65+par->My75*(y[30]+y[31])/NPOP75)
     -y[0]*par->lambdas*(par->Myy*y[4]/NPOPY+par->My25*y[11]/NPOP25+par->My45*y[18]/NPOP45+par->My65*y[25]/NPOP65+par->My75*y[32]/NPOP75);
    f[1] = y[0]*par->lambdaa*(par->Myy*(y[2]+y[3])/NPOPY+par->My25*(y[9]+y[10])/NPOP25+par->My45*(y[16]+y[17])/NPOP45+par->My65*(y[23]+y[24])/NPOP65+par->My75*(y[30]+y[31])/NPOP75)
     +y[0]*par->lambdas*(par->Myy*y[4]/NPOPY+par->My25*y[11]/NPOP25+par->My45*y[18]/NPOP45+par->My65*y[25]/NPOP65+par->My75*y[32]/NPOP75) - par->sigma * y[1];
    f[2] = par->pay *par->sigma * y[1] - par->gammaay * y[2];
    f[3] = (1-par->pay)*par->sigma * y[1] - par->tau * y[3];
    f[4] = par->tau * y[3] - par->deltay * y[4] - par->gammasy * y[4];
    f[5] = par->deltay * y[4] - par->gammaqyused * y[5] - par->ryused * y[5];
    f[6] = par->ryused * y[5];
    // 25-44
    f[7] = -y[7]*par->lambdaa*(par->M25y*(y[2]+y[3])/NPOPY+par->M2525*(y[9]+y[10])/NPOP25+par->M2545*(y[16]+y[17])/NPOP45+par->M2565*(y[23]+y[24])/NPOP65+par->M2575*(y[30]+y[31])/NPOP75)
     -y[7]*par->lambdas*(par->M25y*y[4]/NPOPY+par->M2525*y[11]/NPOP25+par->M2545*y[18]/NPOP45+par->M2565*y[25]/NPOP65+par->M2575*y[32]/NPOP75);
    f[8] = y[7]*par->lambdaa*(par->M25y*(y[2]+y[3])/NPOPY+par->M2525*(y[9]+y[10])/NPOP25+par->M2545*(y[16]+y[17])/NPOP45+par->M2565*(y[23]+y[24])/NPOP65+par->M2575*(y[30]+y[31])/NPOP75)
     +y[7]*par->lambdas*(par->M25y*y[4]/NPOPY+par->M2525*y[11]/NPOP25+par->M2545*y[18]/NPOP45+par->M2565*y[25]/NPOP65+par->M2575*y[32]/NPOP75) - par->sigma * y[8];
    f[9] = par->pa25 *par->sigma * y[8] - par->gammaa25 * y[9];
    f[10] = (1-par->pa25)*par->sigma * y[8] - par->tau * y[10];
    f[11] = par->tau * y[10] - par->delta25 * y[11] - par->gammas25 * y[11];
    f[12] = par->delta25 * y[11] - par->gammaq25used * y[12] - par->r25used * y[12];
    f[13] = par->r25used * y[12];  
    // 45-64
    f[14] = -y[14]*par->lambdaa*(par->M45y*(y[2]+y[3])/NPOPY+par->M4525*(y[9]+y[10])/NPOP25+par->M4545*(y[16]+y[17])/NPOP45+par->M4565*(y[23]+y[24])/NPOP65+par->M4575*(y[30]+y[31])/NPOP75)
     -y[14]*par->lambdas*(par->M45y*y[4]/NPOPY+par->M4525*y[11]/NPOP25+par->M4545*y[18]/NPOP45+par->M4565*y[25]/NPOP65+par->M4575*y[32]/NPOP75);
    f[15] = y[14]*par->lambdaa*(par->M45y*(y[2]+y[3])/NPOPY+par->M4525*(y[9]+y[10])/NPOP25+par->M4545*(y[16]+y[17])/NPOP45+par->M4565*(y[23]+y[24])/NPOP65+par->M4575*(y[30]+y[31])/NPOP75)
     +y[14]*par->lambdas*(par->M45y*y[4]/NPOPY+par->M4525*y[11]/NPOP25+par->M4545*y[18]/NPOP45+par->M4565*y[25]/NPOP65+par->M4575*y[32]/NPOP75) - par->sigma * y[15];
    f[16] = par->pa45 *par->sigma * y[15] - par->gammaa45 * y[16];
    f[17] = (1-par->pa45)*par->sigma * y[15] - par->tau * y[17];
    f[18] = par->tau * y[17] - par->delta45 * y[18] - par->gammas45 * y[18];
    f[19] = par->delta45 * y[18] - par->gammaq45used * y[19] - par->r45used * y[19];
    f[20] = par->r45used * y[19]; 
    // 65-74
    f[21] = -y[21]*par->lambdaa*(par->M65y*(y[2]+y[3])/NPOPY+par->M6525*(y[9]+y[10])/NPOP25+par->M6545*(y[16]+y[17])/NPOP45+par->M6565*(y[23]+y[24])/NPOP65+par->M6575*(y[30]+y[31])/NPOP75)
     -y[21]*par->lambdas*(par->M65y*y[4]/NPOPY+par->M6525*y[11]/NPOP25+par->M6545*y[18]/NPOP45+par->M6565*y[25]/NPOP65+par->M6575*y[32]/NPOP75);
    f[22] = y[21]*par->lambdaa*(par->M65y*(y[2]+y[3])/NPOPY+par->M6525*(y[9]+y[10])/NPOP25+par->M6545*(y[16]+y[17])/NPOP45+par->M6565*(y[23]+y[24])/NPOP65+par->M6575*(y[30]+y[31])/NPOP75)
     +y[21]*par->lambdas*(par->M65y*y[4]/NPOPY+par->M6525*y[11]/NPOP25+par->M6545*y[18]/NPOP45+par->M6565*y[25]/NPOP65+par->M6575*y[32]/NPOP75) - par->sigma * y[22];
    f[23] = par->pa65 *par->sigma * y[22] - par->gammaa65 * y[23];
    f[24] = (1-par->pa65)*par->sigma * y[22] - par->tau * y[24];
    f[25] = par->tau * y[24] - par->delta65 * y[25] - par->gammas65 * y[25];
    f[26] = par->delta65 * y[25] - par->gammaq65used * y[26] - par->r65used * y[26];
    f[27] = par->r65used * y[26];  
    // 75+
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

// function for differential equations nursing homes
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

// social contact data
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

// next generation principle
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

// function for travellers reimportation 
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

// function performing the simulation between specific days (to be called by scenario functions in scenario.c)
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


// initialisation simulation
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

// calling simulation and likelihood 
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

// reading of the data
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


// function for admissible zone for parameters
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




