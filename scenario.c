/*
Corona seiirqd covid-19 program
File: scenario.c (with scenario for simulation)
Copyright Nicolas Franco - UNamur 2020-2021
Public version 1.0 (28/04/2021)
Version corresponding to the scientific paper:
"Covid-19 Belgium: Extended SEIR-QD model with nursing homes and long-term scenarios-based forecasts"
https://doi.org/10.1101/2020.09.07.20190108
*/

/* *************************************************
Fonctions containing the scenarios: to be modified accordingly
-> simucontinu: used for calibration and current behaviour scenario
-> simucontinu2: used for scenario variations (with potential parameters change and change2)
*************************************************
Each period of the simulation must be done with the call of two functions:
->      contact_matrices(par, work, school, leisure,home);
where work, school, leisure,home are proportion of transmission during the period (in [0,1])
->   perform_simu(start,end,home_open,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
where start-end and starting and ending date (day 1 = Nov 29 2020) and home_open=1 if visits are allowed in nursing homes
************************************************* */

#include "corona.h"

void simucontinu(SEIIRQD * etat, SEIIRQD * etat_incid, parameters * par, int end, gsl_rng * r, double reimp[4][DAY_MAX + 1])
{

    int  nhome;
    double  removedfromgenpop = 0;
    double y[40];
    double z[NHOMES_SIMU][NHOMESVAR];
    
    // init

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
// lift of lockdown
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
        if(end<=399){
        contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(360,end,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
        } else {
    contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(360,399,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
     contact_matrices(par, par->work_unlock, 0,par->leisure_sept,par->home_unlock);
    perform_simu(399,415,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
        contact_matrices(par, par->work_unlock, par->school_sept,par->leisure_sept,par->home_unlock);
    perform_simu(415,end,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    
    }

       
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
// lift of lockdown
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
       
// new parameters (to create different scenarios)
  double newwork = par->work_unlock;
  double newleisure = par->leisure_sept;
  double newhome = par->home_unlock;
  newwork = par->work_lock + (par->work_unlock - par->work_lock)*2/4;
  newleisure = par->leisure_lock+(par->leisure_sept-par->leisure_lock)*2/4;
  newhome = par->home_lock + (par->home_unlock - par->home_lock)*2/4;

  double longschool =  par->school_sept; 
  double longwork = par->work_unlock;
  double longleisure = par->leisure_sept;
  double longhome = par->home_unlock;
  longwork = par->work_lock + (par->work_unlock - par->work_lock)*0/4;
  longleisure = par->leisure_lock+(par->leisure_sept-par->leisure_lock)*0/4;
  longhome = par->home_lock + (par->home_unlock - par->home_lock)*0/4;
  longschool = 5.0/6.0*par->school_sept;
  
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
    
            if(end<=399){
    contact_matrices(par, scwork, scschool,scleisure,schome);
    perform_simu(360,end,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
        } else {
    contact_matrices(par, scwork, scschool,scleisure,schome);
    perform_simu(360,399,0,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
           contact_matrices(par, scwork, 0,scleisure,schome);
    perform_simu(399,415,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
        contact_matrices(par, scwork, scschool,scleisure,schome);
    perform_simu(415,end,1,&removedfromgenpop,etat,etat_incid,par,r,y,z,d,dhome,reimp);
    
    }


       
    gsl_odeiv2_driver_free(d);
    gsl_odeiv2_driver_free(dhome);  
}


