//
//  initcomp.h
//  Gibbs
//
//  Created by Yoshi Miyazaki on 2015/02/17.
//  Copyright (c) 2015 Yoshi. All rights reserved.
//

#ifndef Gibbs_initcomp_h
#define Gibbs_initcomp_h

#include "gibbse.h"

using namespace std;

/*--------------------------------------------------------------------------
 Def of class: initcomp
 [What this class does]
 ... obatin initial vector for SA that satisfies
    1. positive mole#
    2. initial composition
     This class basically imports the initial composition data and
     matches with the final composition data. See the example below.
 
 [What to prepare]
 ... input file and considering molecule list data
     Remember to include the molecule in initial composition (input file)
     in the considering molecule list data

 // input file format:
    string(moleculename)     (tab)   double(mole#)
    ...
    ex.
    H                                100.0
    Al                               20.0
    CH30H                            1.5
    O2                               3.0
    SiO2                             0.3
                                       <--- initial_mol (size: # of initial composition mol)
 
 // output : molecule_comp [Vector1d<double> (# of final molecule)]
    H       0.0
    H2     50.0   <---- If input was in element, copy the value to the SIMPLE substance of that element
    Al     20.0
    O2      3.0
    CH3OH   1.5   <---- If input was in molecule, copy the value to the row corresponding the same molecule
    SiO2    0.3
              <--- output_mol (size: # of considering molecule list composition mol)
 ---------------------------------------------------------------------------*/

class initcomp{
public:
    initcomp(string, MoleculeData_G&);
    const Vector1d<double> getmole()   const{return output_mol;};
    const Vector1d<double> getelementcomp() const {return element_mol;};
          tensor1d<double> getmole_t()    {return output_mol.to_tensor();};
    const string           getlist() const{return list_species;};
private:
    MoleculeData     read_initcomp(string);
    
    int              numofcomposition;      // # of composition considering
    string           list_species;          /* list of species */
    Vector1d<double> initial_mol;           // molecule's mol# vector
    Vector1d<double> output_mol;            // element's  mol# vector
    
    Vector1d<double> element_mol;           // Just for now. Not necessary.
};

#endif
