/*
  const.h
  - all the physical constants in the universe

  Created by Yoshi on 2015 Nov. 17
 */

#ifndef astro_const
#define astro_const

/* physical constant */
const double G  = 6.6738e-11;    /* gravitational constant  */
const double Ms = 1.9886e30;     /* solar mass [kg]         */
const double ME = 5.972e24;      /* Earth mass [kg]         */
const double rE = 6371e3;        /* Earth radius [m]        */
const double AU = 1.4959787e11;  /* astronomical unit [m]   */
const double kB = 1.380649e-23;  /* Boltzmann constant      */
const double R  = 8.3144598;     /* gas constant            */
const double mu = 2.2;           /* mean molar mass in disk */
const double mH = 1.673534e-27;  /* atomic hydrogen mass    */
const double NA = 6.02214e23;    /* Avogadro constant       */
const double sb = 5.6704e-8;     /* Stefan-Boltzmann const  */
const double yr2sec = 3.1536e7;
const double Pa2bar = 1e-5;

const double pi = 3.14159265;    /* pi */

/* reference state */
const double Tref = 298.15;     /* reference temperature */
const double Patm = 1.01325e5;  /* reference pressure    */

#endif
