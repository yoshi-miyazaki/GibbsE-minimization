/*
  gibbsminCG.h
  Gibbs energy minimization by conjugate gradient method
  
  by Yoshi Miyazaki on August 4, 2015
  modified on April 8, 2017
  All rights reserved.
*/

#ifndef Gibbs_ConjugateG_h
#define Gibbs_ConjugateG_h

#include <cstdlib>
#include <time.h>
#include <ctime>
#include <float.h>
#include "./const.h"
#include "./tensor.h"
#include "initcomp.h"
#include "solution.h"
#include "massb.h"

const double LIM  = 2.22044604925031e-16*8;
const int print = 0;

/*--------------------------------------------------------------------------
 // Def of class: gibbsminCG
 ... This class provides the equilibrium composition by Gibbs energy minimization,
 .   based on 
 .   - Gibbs energy data `MoleculeData_G` and 
 .   - initial composition defined by `tensor1d<double>`
 .
 .   To use,, 
 .   1. constrcut database:
 .     -> MoleculeData_G  _m(filename, T, P);
 .   2. create gibbsminCG object
 .     -> gibbsminCG      _min(_m, ninit);
 .   3. get result
 .     -> nbest = _min.getnbest();
 .
 .   Functions:
 .   - getnbest():    return the equilibrium composition
 .   - getphase():    return the phase of each species
 .   - convergence(): check the convergence of calculation
 ---------------------------------------------------------------------------*/

class gibbsminCG{
 public:
    gibbsminCG(MoleculeData_G&, tensor1d<double>&);

    /* return optimized values */
    bool             convergence(){return converge;};
    double           getGbest(){return Gbest;};
    tensor1d<int>    getphase(){return phase;};
    tensor1d<double> getnbest(){return nbest;};
    tensor1d<double> getngas_or_solid(int);
    
    /* output result to file */
    void             result(double, double, string&);
    double           meltfrac(double, double, string&);
    double           meltfrac();
    
 private:
    /* pressure and temperature  ... read from MoleculeData_G */
    double           T;
    double           P;
    
    /* calculate Gibbs energy and gradient */
    void             conjugateG();      /* minimize function gibbsE from n_in */
    double           gibbsE(tensor1d<double>&);
    tensor1d<double> grad(tensor1d<double>&, int);
    
    /* solid solution and melt mixing */
    melt_mixing        meltSi;
    tensor1d<solution> ss_list;   /* store the list of solid solution */
    Matrix<double>     ss_matrix;
    void               create_sslist();
    
    /* initial n & massbalance */
    tensor1d<double> n_in;
    tensor1d<double> atomic;
    Matrix<double>   massm;
    
    /* gibbs energy data & phase data */
    MoleculeData_G   list;
    tensor1d<double> gibbse;
    tensor1d<int>    phase;
    tensor1d<double> getGibbslist(MoleculeData_G&);
    tensor1d<int>    getphaselist(MoleculeData_G&);
    
    /* conjugate gradient */
    double           Gbest;
    tensor1d<double> nbest;
    bool             converge;
    
    void             mass_balance(tensor1d<double>&, tensor1d<int>&);
    Matrix<double>   compose_projection(Matrix<double>&);
    void             zerogas(tensor1d<double>&, tensor1d<double>&, tensor1d<int>&);
    void             avoidneg(tensor1d<double>&, tensor1d<double>&, tensor1d<int>&);
    tensor1d<double> shrink  (tensor1d<double>&, tensor1d<double>&, double&);
    tensor1d<double> perturb_zero(tensor1d<double>&, tensor1d<double>&);
    void             bisection(tensor1d<double>, tensor1d<double>&, tensor1d<double>&, double, tensor1d<int>&, int&, int&);
    bool             ck_smaller(tensor1d<double>&, tensor1d<double>&, tensor1d<double>&, tensor1d<int>&);
    double           dG_direction(tensor1d<double>&, tensor1d<double>);

    /* projection matrix P = I - B (BT B)-1 BT */
    Matrix<double>   projv;
    Matrix<double>   projv_saved;    /* reduce calculation */
    tensor1d<int>    exist_saved;    /* if diminish is same as previous */
};

#endif
