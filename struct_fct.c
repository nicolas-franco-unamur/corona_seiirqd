/*
Corona seiirqd covid-19 program
File: struct_fct.c (to get element from struc defined in corona.h)
Copyright Nicolas Franco - UNamur 2020-2021
Public version 1.0 (28/04/2021)
Version corresponding to the scientific paper:
"Covid-19 Belgium: Extended SEIR-QD model with nursing homes and long-term scenarios-based forecasts"
https://doi.org/10.1101/2020.09.07.20190108
*/

/*
Functions catching element of the different stuctures: should be modified if structure in corona.h are modified
*/

#include "corona.h"

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

// basic function for seed
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
