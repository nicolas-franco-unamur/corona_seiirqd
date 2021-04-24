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
    double S;			//2 - 65
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
    double p0;			//1
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


// loop through elements
double *foreachetat(int i, SEIIRQD* etat){
switch(i){
default:
exit(-23);
case 1: return &etat->S;
case 2: return &etat->E;
case 3: return &etat->AI;
case 4: return &etat->PI;
case 5: return &etat->SI;
case 6: return &etat->R;
case 7: return &etat->Q;
case 8: return &etat->D;
case 9: return &etat->SY;
case 10: return &etat->EY;
case 11: return &etat->AIY;
case 12: return &etat->PIY;
case 13: return &etat->SIY;
case 14: return &etat->RY;
case 15: return &etat->QY;
case 16: return &etat->DY;
case 17: return &etat->S25;
case 18: return &etat->E25;
case 19: return &etat->AI25;
case 20: return &etat->PI25;
case 21: return &etat->SI25;
case 22: return &etat->R25;
case 23: return &etat->Q25;
case 24: return &etat->D25;
case 25: return &etat->S45;
case 26: return &etat->E45;
case 27: return &etat->AI45;
case 28: return &etat->PI45;
case 29: return &etat->SI45;
case 30: return &etat->R45;
case 31: return &etat->Q45;
case 32: return &etat->D45;
case 33: return &etat->S65;
case 34: return &etat->E65;
case 35: return &etat->AI65;
case 36: return &etat->PI65;
case 37: return &etat->SI65;
case 38: return &etat->R65;
case 39: return &etat->Q65;
case 40: return &etat->D65;
case 41: return &etat->S75;
case 42: return &etat->E75;
case 43: return &etat->AI75;
case 44: return &etat->PI75;
case 45: return &etat->SI75;
case 46: return &etat->R75;
case 47: return &etat->Q75;
case 48: return &etat->D75;
case 49: return &etat->SH;
case 50: return &etat->EH;
case 51: return &etat->AIH;
case 52: return &etat->PIH;
case 53: return &etat->SIH;
case 54: return &etat->RH;
case 55: return &etat->QH;
case 56: return &etat->DH;
case 57: return &etat->D75H;
case 58: return &etat->SIQ;
case 59: return &etat->HQ;
case 60: return &etat->QR;
case 61: return &etat->QD;
case 62: return &etat->AIR;
case 63: return &etat->SIR;
}
}

double *foreachetatarray(int i, SEIIRQDarray* etat){
switch(i){
default:
exit(-23);
case 1: return etat->S;
case 2: return etat->E;
case 3: return etat->AI;
case 4: return etat->PI;
case 5: return etat->SI;
case 6: return etat->R;
case 7: return etat->Q;
case 8: return etat->D;
case 9: return etat->SY;
case 10: return etat->EY;
case 11: return etat->AIY;
case 12: return etat->PIY;
case 13: return etat->SIY;
case 14: return etat->RY;
case 15: return etat->QY;
case 16: return etat->DY;
case 17: return etat->S25;
case 18: return etat->E25;
case 19: return etat->AI25;
case 20: return etat->PI25;
case 21: return etat->SI25;
case 22: return etat->R25;
case 23: return etat->Q25;
case 24: return etat->D25;
case 25: return etat->S45;
case 26: return etat->E45;
case 27: return etat->AI45;
case 28: return etat->PI45;
case 29: return etat->SI45;
case 30: return etat->R45;
case 31: return etat->Q45;
case 32: return etat->D45;
case 33: return etat->S65;
case 34: return etat->E65;
case 35: return etat->AI65;
case 36: return etat->PI65;
case 37: return etat->SI65;
case 38: return etat->R65;
case 39: return etat->Q65;
case 40: return etat->D65;
case 41: return etat->S75;
case 42: return etat->E75;
case 43: return etat->AI75;
case 44: return etat->PI75;
case 45: return etat->SI75;
case 46: return etat->R75;
case 47: return etat->Q75;
case 48: return etat->D75;
case 49: return etat->SH;
case 50: return etat->EH;
case 51: return etat->AIH;
case 52: return etat->PIH;
case 53: return etat->SIH;
case 54: return etat->RH;
case 55: return etat->QH;
case 56: return etat->DH;
case 57: return etat->D75H;
case 58: return etat->SIQ;
case 59: return etat->HQ;
case 60: return etat->QR;
case 61: return etat->QD;
case 62: return etat->AIR;
case 63: return etat->SIR;

}
}



double *foreachpar(int i,parameters* par){
switch(i){
default:
exit(-23);
case 1: return &par->p0;
case 2: return &par->lambdaa;
case 3: return &par->lambdas;
case 4: return &par-> sigma;
case 5: return &par-> tau;
case 6: return &par-> pay;
case 7: return &par-> pa25;
case 8: return &par-> pa45;
case 9: return &par-> pa65;
case 10: return &par-> pa75;
case 11: return &par-> pah;
case 12: return &par-> deltay;
case 13: return &par-> delta25;
case 14: return &par-> delta45;
case 15: return &par-> delta65;
case 16: return &par-> delta75;
case 17: return &par-> deltah;
case 18: return &par-> gammaay;
case 19: return &par-> gammaa25;
case 20: return &par-> gammaa45;
case 21: return &par-> gammaa65;
case 22: return &par-> gammaa75;
case 23: return &par-> gammaah;
case 24: return &par-> gammasy;
case 25: return &par-> gammas25;
case 26: return &par-> gammas45;
case 27: return &par-> gammas65;
case 28: return &par-> gammas75;
case 29: return &par-> gammash;
case 30: return &par-> gammaqy;
case 31: return &par-> gammaq25;
case 32: return &par-> gammaq45;
case 33: return &par-> gammaq65;
case 34: return &par-> gammaq75;
case 35: return &par-> gammaqh;
case 36: return &par-> ry;
case 37: return &par-> r25;
case 38: return &par-> r45;
case 39: return &par-> r65;
case 40: return &par-> r75;
case 41: return &par-> rh;
case 42: return &par-> rht;
case 43: return &par-> recovery;
case 44: return &par-> rcap;
case 45: return &par-> rcaps;
case 46: return &par-> supphosp;
case 47: return &par-> cap;
case 48: return &par-> caps;
case 49: return &par-> delay;
case 50: return &par-> pconf;
case 51: return &par-> pth;
case 52: return &par-> pthp;
case 53: return &par-> mh;
case 54: return &par-> home_lock;
case 55: return &par-> work_lock;
case 56: return &par-> leisure_lock;
case 57: return &par-> home_unlock;
case 58: return &par-> work_unlock;
case 59: return &par-> school_unlock;
case 60: return &par-> leisure_june;
case 61: return &par-> leisure_july;
case 62: return &par-> leisure_augustus;
case 63: return &par-> school_sept;
case 64: return &par-> reimp;
case 65: return &par-> leisure_sept;
}
}


double *foreachpararray(int i,parametersarray* par){
switch(i){
default:
exit(-23);
case 1: return par->p0;
case 2: return par->lambdaa;
case 3: return par->lambdas;
case 4: return par-> sigma;
case 5: return par-> tau;
case 6: return par-> pay;
case 7: return par-> pa25;
case 8: return par-> pa45;
case 9: return par-> pa65;
case 10: return par-> pa75;
case 11: return par-> pah;
case 12: return par-> deltay;
case 13: return par-> delta25;
case 14: return par-> delta45;
case 15: return par-> delta65;
case 16: return par-> delta75;
case 17: return par-> deltah;
case 18: return par-> gammaay;
case 19: return par-> gammaa25;
case 20: return par-> gammaa45;
case 21: return par-> gammaa65;
case 22: return par-> gammaa75;
case 23: return par-> gammaah;
case 24: return par-> gammasy;
case 25: return par-> gammas25;
case 26: return par-> gammas45;
case 27: return par-> gammas65;
case 28: return par-> gammas75;
case 29: return par-> gammash;
case 30: return par-> gammaqy;
case 31: return par-> gammaq25;
case 32: return par-> gammaq45;
case 33: return par-> gammaq65;
case 34: return par-> gammaq75;
case 35: return par-> gammaqh;
case 36: return par-> ry;
case 37: return par-> r25;
case 38: return par-> r45;
case 39: return par-> r65;
case 40: return par-> r75;
case 41: return par-> rh;
case 42: return par-> rht;
case 43: return par-> recovery;
case 44: return par-> rcap;
case 45: return par-> rcaps;
case 46: return par-> supphosp;
case 47: return par-> cap;
case 48: return par-> caps;
case 49: return par-> delay;
case 50: return par-> pconf;
case 51: return par-> pth;
case 52: return par-> pthp;
case 53: return par-> mh;
case 54: return par-> home_lock;
case 55: return par-> work_lock;
case 56: return par-> leisure_lock;
case 57: return par-> home_unlock;
case 58: return par-> work_unlock;
case 59: return par-> school_unlock;
case 60: return par-> leisure_june;
case 61: return par-> leisure_july;
case 62: return par-> leisure_augustus;
case 63: return par-> school_sept;
case 64: return par-> reimp;
case 65: return par-> leisure_sept;
}
}

const char* foreachparname(int i){
switch(i){
default:
exit(-23);
case 1: return "p0";
case 2: return "lambdaa";
case 3: return "lambdas";
case 4: return "sigma";
case 5: return "tau";
case 6: return "pay";
case 7: return "pa25";
case 8: return "pa45";
case 9: return "pa65";
case 10: return "pa75";
case 11: return "pah";
case 12: return "deltay";
case 13: return "delta25";
case 14: return "delta45";
case 15: return "delta65";
case 16: return "delta75";
case 17: return "deltah";
case 18: return "gammaay";
case 19: return "gammaa25";
case 20: return "gammaa45";
case 21: return "gammaa65";
case 22: return "gammaa75";
case 23: return "gammaah";
case 24: return "gammasy";
case 25: return "gammas25";
case 26: return "gammas45";
case 27: return "gammas65";
case 28: return "gammas75";
case 29: return "gammash";
case 30: return "gammaqy";
case 31: return "gammaq25";
case 32: return "gammaq45";
case 33: return "gammaq65";
case 34: return "gammaq75";
case 35: return "gammaqh";
case 36: return "ry";
case 37: return "r25";
case 38: return "r45";
case 39: return "r65";
case 40: return "r75";
case 41: return "rh";
case 42: return "rht";
case 43: return "recovery";
case 44: return "rcap";
case 45: return "rcaps";
case 46: return "supphosp";
case 47: return "cap";
case 48: return "caps";
case 49: return "delay";
case 50: return "pconf";
case 51: return "pth";
case 52: return "pthp";
case 53: return "mh";
case 54: return "home_lock";
case 55: return "work_lock";
case 56: return "leisure_lock";
case 57: return "home_unlock";
case 58: return "work_unlock";
case 59: return "school_unlock";
case 60: return "leisure_june";
case 61: return "leisure_july";
case 62: return "leisure_augustus";
case 63: return "school_sept";
case 64: return "reimp";
case 65: return "leisure_sept";
}
}

