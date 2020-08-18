//  massb_oxide.h
//
//  Created by Yoshi Miyazaki on 2017 Mar 7.
//  Copyright (c) 2017 Yoshi. All rights reserved.

#ifndef Gibbs_massBalance_h
#define Gibbs_massBalance_h

#include "gibbse.h"

const int is_oxide = false;

// for minimum norm & null space
class massBalance{
public:
    massBalance(MoleculeData_G& moleculelist);
    
    Matrix<double>   get_massBalance(){ return massm; };
    tensor1d<double> get_atomic()     { return atomic;};
    
private:
    int              d;        /* # of elements in the system  */
    int              m;        /* # of species considered here */
    Matrix<double>   massm;    /* mass-balance matrix          */
    tensor1d<double> atomic;   /* represents element # of each column */
    
    Matrix<double>   nspace;
    tensor1d<int>    which;  // indicate relation nspace & molecule

    void   create_massBalance(MoleculeData_G&);
    void   create_nullspace(MoleculeData_G&);

    int    getdim(MoleculeData_G&);
    int    issimple(Molecule&);
    int    simple_index(const int);
    string simple(const int);
};


#endif
