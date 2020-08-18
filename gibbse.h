//
//  gibbse.h
//  Gibbs
//
//  Created by Yoshi Miyazaki on 2015/01/06.
//  Copyright (c) 2015 Yoshi. All rights reserved.
//

#ifndef Gibbs_gibbse_h
#define Gibbs_gibbse_h

#include <cmath>
#include <sstream>
#include <tuple>
#include <ctime>
#include "element.h"
#include "./const.h"

/*--------------------------------------------------------------------------
 Def of class: Molecule_Gibbs
   ---> inheritance class of Molecule
 [What this class does]
 ... obatin apparent Gibbs energy & phase information
     for temperature T and pressure P
 [What to prepare]
 ... thermodynamics table (extracted from JANAF)
     use converted JANAF table using JANAF.cpp
 T       CP     S       GHr/T   HHr/T       ∆fH     ∆fG     logKf
 298     8.899	-94.05	51.066	-109.275	80.112  80.112  10.00
 ---------------------------------------------------------------------------*/

class Molecule_Gibbs : public Molecule{
public:
    Molecule_Gibbs(string, double, double);
    Molecule_Gibbs(string, double, double, Vector1d<double>&);
    double getGibbsE(){return gibbsE;};
    int    getphase() {return phase;};
    
    double getT()  { return T;};
    double getP()  { return P;};
    double getK0() { return K0;};
    double getK0d(){ return K0d;};
private:
    double gibbsE;
    double K0;
    double K0d;
    double T;
    double P;
    int    phase;                                    // gas=0, liquid=1, solid=2
    
    double data_search(string);
    double GibbsE_Stixrude(double, double, double, double, double, double, double);
    double GibbsE_EoS(double, double, double, double, double, double, double, double, double);
    double GibbsE_melting(ifstream&);
    tuple<double, char> GibbsE_JANAF(ifstream&, string);
    
    /* basic functions */
    void   set_phase(char);
    double internal_division(double, double, double, double, double);
    
    /* functions for Stixrude */
    double GibbsE_solid(string, double, double);
    double debyeF_func(double, double, double, double, double);
    
    /* fucntions for melt*/
    double int_VdP(double, double, double, double, double);
    double int_dSdP(double, double, double, double, double, double);
    double dGdT(string, double, double);
    double dGdP(string, double, double);

    /* functions not used */
    double dH1(double, double, double, double, double, double);
    double dS1(double, double, double, double, double, double);
};

/*--------------------------------------------------------------------------
 Def of class: MoleculeData_Gibbs
   ---> contaning Molecule_Gibbs information
 [What this class does]
 ... Read file contaning molecule names and make Molecule_Gibbs object
     for each name
 [What to prepare]
 ... file.txt conataning the list of molecules
     ex.) H2, H2O, SiO2+4H2O, Mg1.8Fe0.2SiO4
 ---------------------------------------------------------------------------*/
class MoleculeData_G{
public:
    MoleculeData_G(){};
    MoleculeData_G(string, double, double);
    Molecule_Gibbs& operator[](const int i){return database[i];};  /* Return i'th molecule */
    
    int    intspec(string);
    int    getnumofMolecule(){ return int(numofMolecule);};
    string getlist()               { return list_species; };
    Vector1d<double> getmolarmass(){ return molarmass;    };
    tensor1d<double> getmolarmass_t(){ return molarmass.to_tensor(); };

    MoleculeData_G& operator=(const MoleculeData_G&);  /* assign */
    
    void rm_nonexisting(Vector1d<double>&);            /* delete species that includes elements not used */
    void rm_nonexisting(tensor1d<double>&);            /* delete species that includes elements not used */
    
private:
    size_t numofMolecule;
    string list_species;
    vector<Molecule_Gibbs> database;
    Vector1d<double>       molarmass;
    
    void   erase(int);
};


#endif
