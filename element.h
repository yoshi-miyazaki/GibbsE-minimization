//
//  element.h
//  Gibbs
//
//  Created by Yoshi Miyazaki on 2014/10/15.
//  Copyright (c) 2014 Yoshi. All rights reserved.
//

#ifndef Gibbs_element_h
#define Gibbs_element_h

#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
#include "./matrix.h"
using namespace std;

const int numofElement = 30;

/*--------------------------------------------------------------------------
 // Def of class: Element
 ---------------------------------------------------------------------------*/
class Element{
public:
    Element(){};
    void setNumber(int _i){AtomicNumber= _i;};
    int  getNumber(){return AtomicNumber;};
    
    void setName(string _s){ElementName= _s;};
    void showElementName(){cout << "Element is " << ElementName << endl;};
    string getElementName(){return ElementName;}
    bool isthisElement(string _s){return _s==ElementName ? 1 : 0;};
    
    void setWeight(double _i){AtomicWeight= _i;};
    double getWeight(){return AtomicWeight;}
private:
    int    AtomicNumber;
    string ElementName;
    double AtomicWeight;
};

/*--------------------------------------------------------------------------
 // Def of class: ElementData
 ---------------------------------------------------------------------------*/
class ElementData{
private:
    size_t totalNumofElement;
    vector<Element> ElementDatabase;
public:
    ElementData();                               // Database written in the cpp file
    ~ElementData();
    Element operator[](const int i);             // Return i'th element
    Element correspond(string _s);
    size_t getElementDataSize(){return totalNumofElement;};
};

/*--------------------------------------------------------------------------
 // Def of class: Molecule
 ---------------------------------------------------------------------------*/
class Molecule{
public:
    Molecule(){};
    Molecule(string _s);
    void   setName(string _s){MoleculeName= _s;};
    void   showMoleculeName();
    Vector1d<double> getComposition(){return composition;};
    double getComposition(int i){return composition[i];};
    double getWeight()      {return MolecularWeight;};
    int    getno_atoms()    {return no_atoms;}
    string getMoleculeName(){return MoleculeName;};
    Molecule& operator=(const Molecule &assign);     // Assignment
    bool   operator==(const Molecule &same)
                {if (MoleculeName==same.MoleculeName){return 1;} else{return 0;}};       // same
protected:
    string MoleculeName;
    Vector1d<double> composition;
    int    no_atoms;
    double MolecularWeight;
};

/*--------------------------------------------------------------------------
 Def of class: MoleculeData
 ---------------------------------------------------------------------------*/
class MoleculeData{
public:
    MoleculeData(){};
    MoleculeData(string);
    Molecule& operator[](const int i){return MoleculeDatabase[i];};             // Return i'th molecule
    int    getnumofElement(){return int(numofElement);};
    int    getnumofMolecule(){return int(numofMolecule);};
    string getlist(){return list_species;}
private:
    size_t numofElement;
    size_t numofMolecule;
    string list_species;
    vector<Molecule> MoleculeDatabase;
};

#endif
