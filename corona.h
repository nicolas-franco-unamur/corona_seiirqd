/*
Corona seiirqd covid-19 program
File: corona.h (header)
Copyright Nicolas Franco - UNamur 2020-2021
Public version 1.0 (28/04/2021)
Version corresponding to the scientific paper:
"Covid-19 Belgium: Extended SEIR-QD model with nursing homes and long-term scenarios-based forecasts"
https://doi.org/10.1101/2020.09.07.20190108
*/

/*
Main constants and definition of the compartiments: should not be modified except for an important modification of the model
*/

#ifndef CORONA_H_
#define CORONA_H_

#include "command.h"

// librairies

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

// constants

#define NPOP    11500000
#define NPOPY   3250000
#define NPOP25  3000000
#define NPOP45  3080000
#define NPOP65  1150000
#define NPOP75H 1020000
#define NPOPH   150000
#define NHOMES    (2000)
#define NPOPPERH  (NPOPH/NHOMES)
#define NPOP75 (NPOP75H-NPOPH)
#if SPEED > 1
#define NHOMES_SIMU (NHOMES/SPEED)
#define COEFSPEED SPEED
#else
#define NHOMES_SIMU NHOMES
#define COEFSPEED 1
#endif

#define NHOMESVAR 11
#define NUMSEG 63
#define NUMPAR 65

// prior + sigma gaussian
#define PO 0.0002
#define PO_SGM   0.0000002
#define LAMBDA 0.08
#define LAMBDA_SGM 0.00005
#define SIGMA 0.5
#define SIGMA_SGM 0.0005
#define TAU 0.2
#define TAU_SGM 0.0002
#define PA 0.6
#define PA_SGM  0.0005
#define PA_DELTA 0.1
#define DELTA 0.05
#define DELTA_SGM  0.00005
#define DELTA_DELTA 0.005
#define GAMMA 0.25
#define GAMMA_SGM  0.0002
#define GAMMA_DELTA 0.02
#define GAMMAQ 0.05
#define GAMMAQ_SGM  0.00005
#define GAMMAQ_DELTA 0.01
#define RR 0.02
#define RR_SGM  0.00002
#define RR_DELTA 0.005
#define CAP 4000
#define CAP_SGM  5
#define CAPS 2000
#define CAPS_SGM  2
#define DELAY 15
#define DELAY_SGM 0.02
#define PCONF 0.8
#define PCONF_SGM 0.0005
#define PTH 10.0
#define PTH_SGM 0.01
#define PTHP 10.0
#define PTHP_SGM 0.01
#define MH 0.5
#define MH_SGM  0.0005
#define SUPPHOSP 1.15
#define CONTACT 0.1
#define CONTACT_HOME 0.5
#define CONTACT_SGM 0.0001
#define NEWCONTACT_SGM 0.0005
#define CONTACT_DELTA 0.05
#define RECOVERY 0.7
#define RECOVERY_SGM 0.0005
#define REIMP 50
#define REIMP_SGM 0.05

typedef struct SEIIRQD {
    double S;			//2 - 65  // number of the compartment in the output file
    double E;			//3 - 66
    double AI;			//4 - 67
    double PI;			//5 - 68
    double SI;			//6 - 69
    double R;			//7 - 70
    double Q;			//8 - 71
    double D;			//9 - 72
    double SY;			//10 - 
    double EY;			//11 - 
    double AIY;			//12 - 
    double PIY;			//13 - 
    double SIY;			//14 - 
    double RY;			//15 - 
    double QY;			//16 - 
    double DY;			//17 - 80
    double S25;			//18 - 
    double E25;			//19 - 
    double AI25;	    //20 - 
    double PI25;		//21 - 
    double SI25;		//22 - 
    double R25;			//23 - 
    double Q25;			//24 - 
    double D25;			//25 - 88
    double S45;			//26 - 
    double E45;			//27 - 
    double AI45;		//28 - 
    double PI45;		//29 - 
    double SI45;		//30 - 
    double R45;			//31 - 
    double Q45;			//32 - 
    double D45;			//33 - 96
    double S65;			//34 - 
    double E65;			//35 - 
    double AI65;		//36 - 
    double PI65;		//37 - 
    double SI65;		//38 - 
    double R65;			//39 - 
    double Q65;			//40 - 
    double D65;			//41 - 104
    double S75;			//42 - 
    double E75;			//43 - 
    double AI75;		//44 - 
    double PI75;		//45 - 
    double SI75;		//46 - 
    double R75;			//47 - 
    double Q75;			//48 - 
    double D75;			//49 - 112
    double SH;			//50 - 
    double EH;			//51 - 
    double AIH;			//52 - 
    double PIH;			//53 - 
    double SIH;			//54 - 
    double RH;			//55 - 
    double QH;			//56 - 
    double DH;			//57 - 
    double D75H;		//58 - 121
    double SIQ;			//59 - 122
    double HQ;			//60 - 123
    double QR;			//61 - 124
    double QD;			//62 - 125
    double AIR;			//63 - 126
    double SIR;			//64 - 127
} SEIIRQD;

typedef struct SEIIRQDarray {
    double S[SIZE_SAMPLE];			//2
    double E[SIZE_SAMPLE];			//3
    double AI[SIZE_SAMPLE];			//4
    double PI[SIZE_SAMPLE];			//5
    double SI[SIZE_SAMPLE];			//6
    double R[SIZE_SAMPLE];			//7
    double Q[SIZE_SAMPLE];			//8
    double D[SIZE_SAMPLE];			//9
    double SY[SIZE_SAMPLE];			//10
    double EY[SIZE_SAMPLE];			//11
    double AIY[SIZE_SAMPLE];			//12
    double PIY[SIZE_SAMPLE];			//13
    double SIY[SIZE_SAMPLE];			//14
    double RY[SIZE_SAMPLE];			//15
    double QY[SIZE_SAMPLE];			//16
    double DY[SIZE_SAMPLE];			//17
    double S25[SIZE_SAMPLE];			//18
    double E25[SIZE_SAMPLE];			//19
    double AI25[SIZE_SAMPLE];	    //20
    double PI25[SIZE_SAMPLE];		//21
    double SI25[SIZE_SAMPLE];		//22
    double R25[SIZE_SAMPLE];			//23
    double Q25[SIZE_SAMPLE];			//24
    double D25[SIZE_SAMPLE];			//25
    double S45[SIZE_SAMPLE];			//26
    double E45[SIZE_SAMPLE];			//27
    double AI45[SIZE_SAMPLE];		//28
    double PI45[SIZE_SAMPLE];		//29
    double SI45[SIZE_SAMPLE];		//30
    double R45[SIZE_SAMPLE];			//31
    double Q45[SIZE_SAMPLE];			//32
    double D45[SIZE_SAMPLE];			//33
    double S65[SIZE_SAMPLE];			//34
    double E65[SIZE_SAMPLE];			//35
    double AI65[SIZE_SAMPLE];		//36
    double PI65[SIZE_SAMPLE];		//37
    double SI65[SIZE_SAMPLE];		//38
    double R65[SIZE_SAMPLE];			//39
    double Q65[SIZE_SAMPLE];			//40
    double D65[SIZE_SAMPLE];			//41
    double S75[SIZE_SAMPLE];			//42
    double E75[SIZE_SAMPLE];			//43
    double AI75[SIZE_SAMPLE];		//44
    double PI75[SIZE_SAMPLE];		//45
    double SI75[SIZE_SAMPLE];		//46
    double R75[SIZE_SAMPLE];			//47
    double Q75[SIZE_SAMPLE];			//48
    double D75[SIZE_SAMPLE];			//49
    double SH[SIZE_SAMPLE];			//50
    double EH[SIZE_SAMPLE];			//51
    double AIH[SIZE_SAMPLE];			//52
    double PIH[SIZE_SAMPLE];			//53
    double SIH[SIZE_SAMPLE];			//54
    double RH[SIZE_SAMPLE];			//55
    double QH[SIZE_SAMPLE];			//56
    double DH[SIZE_SAMPLE];			//57
    double D75H[SIZE_SAMPLE];		//58
    double SIQ[SIZE_SAMPLE];			//59
    double HQ[SIZE_SAMPLE];			//60
    double QR[SIZE_SAMPLE];			//61
    double QD[SIZE_SAMPLE];			//62
    double AIR[SIZE_SAMPLE];			//63
    double SIR[SIZE_SAMPLE];			//64
} SEIIRQDarray;


typedef struct parameters {
    double p0;			//1  // number of the parameter in the output file
    double lambdaa;     //2
    double lambdas;     //3
    double sigma;		//4
    double tau;	    	//5
    double pay;         //6
    double pa25;        //7
    double pa45;        //8
    double pa65;        //9
    double pa75;        //10
    double pah;         //11
    double deltay;		//12
    double delta25;		//13
    double delta45;		//14
    double delta65;		//15
    double delta75;		//16
    double deltah;		//17
    double gammaay;		//18
    double gammaa25;	//19
    double gammaa45;	//20
    double gammaa65;	//21
    double gammaa75;	//22
    double gammaah;		//23
    double gammasy;		//24
    double gammas25;	//25
    double gammas45;	//26
    double gammas65;	//27
    double gammas75;	//28
    double gammash;		//29
    double gammaqy;		//30
    double gammaq25;	//31
    double gammaq45;	//32
    double gammaq65;	//33
    double gammaq75;	//34
    double gammaqh;		//35
    double ry;			//36
    double r25;			//37
    double r45;			//38
    double r65;			//39
    double r75;			//40
    double rh;			//41
    double rht;			//42
    double recovery;    //43
    double rcap;        //44
    double rcaps;       //45
    double supphosp;    //46
    double cap;			//47
    double caps;		//48
    double delay;       //49
    double pconf;		//50
    double pth;			//51
    double pthp;		//52
    double mh;		    //53
    double home_lock;   //54
    double work_lock;   // 55
    double leisure_lock;// 56
    double home_unlock; //57
    double work_unlock; // 58
    double school_unlock;    // 59
    double leisure_june;    // 60
    double leisure_july;    // 61
    double leisure_augustus; // 62
    double school_sept; // 63
    double reimp; // 64
    double leisure_sept; // 65
        // additional temp parameters
    double Myy;
    double My25;
    double My45;
    double My65;
    double My75;
    double M25y;
    double M2525;
    double M2545;
    double M2565;
    double M2575;
    double M45y;
    double M4525;
    double M4545;
    double M4565;
    double M4575;
    double M65y;
    double M6525;
    double M6545;
    double M6565;
    double M6575;
    double M75y;
    double M7525;
    double M7545;
    double M7565;
    double M7575;
    double deltahused;
    double rhtused;
    double gammaqyused;		//30
    double gammaq25used;	//31
    double gammaq45used;	//32
    double gammaq65used;	//33
    double gammaq75used;	//34
    double gammaqhused;		//35
    double ryused;			//36
    double r25used;			//37
    double r45used;			//38
    double r65used;			//39
    double r75used;			//40
    double rhused;			//41

} parameters;


typedef struct parametersarray {
    double p0[SIZE_SAMPLE];			//1
    double lambdaa[SIZE_SAMPLE];     //2
    double lambdas[SIZE_SAMPLE];     //3
    double sigma[SIZE_SAMPLE];		//4
    double tau[SIZE_SAMPLE];	    	//5
    double pay[SIZE_SAMPLE];         //6
    double pa25[SIZE_SAMPLE];        //7
    double pa45[SIZE_SAMPLE];        //8
    double pa65[SIZE_SAMPLE];        //9
    double pa75[SIZE_SAMPLE];        //10
    double pah[SIZE_SAMPLE];         //11
    double deltay[SIZE_SAMPLE];		//12
    double delta25[SIZE_SAMPLE];		//13
    double delta45[SIZE_SAMPLE];		//14
    double delta65[SIZE_SAMPLE];		//15
    double delta75[SIZE_SAMPLE];		//16
    double deltah[SIZE_SAMPLE];		//17
    double gammaay[SIZE_SAMPLE];		//18
    double gammaa25[SIZE_SAMPLE];	//19
    double gammaa45[SIZE_SAMPLE];	//20
    double gammaa65[SIZE_SAMPLE];	//21
    double gammaa75[SIZE_SAMPLE];	//22
    double gammaah[SIZE_SAMPLE];		//23
    double gammasy[SIZE_SAMPLE];		//24
    double gammas25[SIZE_SAMPLE];	//25
    double gammas45[SIZE_SAMPLE];	//26
    double gammas65[SIZE_SAMPLE];	//27
    double gammas75[SIZE_SAMPLE];	//28
    double gammash[SIZE_SAMPLE];		//29
    double gammaqy[SIZE_SAMPLE];		//30
    double gammaq25[SIZE_SAMPLE];	//31
    double gammaq45[SIZE_SAMPLE];	//32
    double gammaq65[SIZE_SAMPLE];	//33
    double gammaq75[SIZE_SAMPLE];	//34
    double gammaqh[SIZE_SAMPLE];		//35
    double ry[SIZE_SAMPLE];			//36
    double r25[SIZE_SAMPLE];			//37
    double r45[SIZE_SAMPLE];			//38
    double r65[SIZE_SAMPLE];			//39
    double r75[SIZE_SAMPLE];			//40
    double rh[SIZE_SAMPLE];			//41
    double rht[SIZE_SAMPLE];			//42
    double recovery[SIZE_SAMPLE];    //43
    double rcap[SIZE_SAMPLE];        //44
    double rcaps[SIZE_SAMPLE];       //45
    double supphosp[SIZE_SAMPLE];    //46
    double cap[SIZE_SAMPLE];			//47
    double caps[SIZE_SAMPLE];		//48
    double delay[SIZE_SAMPLE];       //49
    double pconf[SIZE_SAMPLE];		//50
    double pth[SIZE_SAMPLE];			//51
    double pthp[SIZE_SAMPLE];		//52
    double mh[SIZE_SAMPLE];		    //53
    double home_lock[SIZE_SAMPLE];   //54
    double work_lock[SIZE_SAMPLE];   // 55
    double leisure_lock[SIZE_SAMPLE];// 56
    double home_unlock[SIZE_SAMPLE]; //57
    double work_unlock[SIZE_SAMPLE]; // 58
    double school_unlock[SIZE_SAMPLE];    // 59
    double leisure_june[SIZE_SAMPLE];    // 60
    double leisure_july[SIZE_SAMPLE];    // 61
    double leisure_augustus[SIZE_SAMPLE]; // 62
    double school_sept[SIZE_SAMPLE]; // 63
    double reimp[SIZE_SAMPLE]; // 64
    double leisure_sept[SIZE_SAMPLE]; // 65
} parametersarray;


// fct declarations

double *foreachetat(int i, SEIIRQD* etat);
double *foreachetatarray(int i, SEIIRQDarray* etat);
double *foreachpar(int i,parameters* par);
double *foreachpararray(int i,parametersarray* par);
const char* foreachparname(int i);
double likelihood(SEIIRQD * etat, SEIIRQD * etat_incid, SEIIRQD * data, SEIIRQD * data_incid, parameters * par);
int func(double t, const double y[], double f[], void *params);
int funchome(double t, const double y[], double f[], void *params);
void contact_matrices(parameters * par, double work, double school, double leisure, double home);
double compute_rzero(parameters * par, double work, double school, double leisure, double home);
double reimpfct(int day, double reimp[4][DAY_MAX + 1]);
void perform_simu(int start, int end, int openhomes, double *removedfromgenpop,SEIIRQD * etat, SEIIRQD * etat_incid, parameters * par, gsl_rng * r, double *y, double z[NHOMES_SIMU][NHOMESVAR], gsl_odeiv2_driver * d, gsl_odeiv2_driver * dhome, double reimp[4][DAY_MAX + 1]);
void simucontinu(SEIIRQD * etat, SEIIRQD * etat_incid, parameters * par, int end, gsl_rng * r, double reimp[4][DAY_MAX + 1]);
void simucontinu2(SEIIRQD * etat, SEIIRQD * etat_incid, parameters * par, int end, double change,double change2, gsl_rng * r, double reimp[4][DAY_MAX + 1]);
void init(SEIIRQD * etat, SEIIRQD * etat_incid, parameters * par);
double simulikelihood(SEIIRQD * etat, SEIIRQD * etat_incid, SEIIRQD * data, SEIIRQD * data_incid, parameters * par, gsl_rng * r, double reimp[4][DAY_MAX + 1]);
void readdata(SEIIRQD * data, SEIIRQD * data_incid, double reimp[4][DAY_MAX + 1]);
int test_par(parameters par);
parameters random_step_all(gsl_rng * r);
parameters random_step_one(gsl_rng * r);
parameters randomize(parameters temp_par, gsl_rng * r, double coef);
parameters opti_randomize(parameters parburninsamples, SEIIRQD * data, SEIIRQD * data_incid, double *burninlikelihood, int burnin, gsl_rng * r,double reimp[4][DAY_MAX + 1]);
unsigned long int random_seed();
void output( gsl_rng * r,double reimp[4][DAY_MAX + 1]);


#endif