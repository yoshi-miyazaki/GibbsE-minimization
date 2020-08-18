//
//  element.cpp
//  Gibbs
//
//  Created by Yoshi Miyazaki on 2015/04/11.
//  Copyright (c) 2015 Yoshi. All rights reserved.
//

#include "element.h"

/*--------------------------------------------------------------------------
 class: Molecule
 ---------------------------------------------------------------------------*/
Molecule::Molecule(string name){
    MoleculeName = name;
    
    // Element Database
    ElementData ElementDatabase;
    
    // Construct Composition Vector (size : num of Elements considered)
    composition.resize(0.0, numofElement);
    
    // Molecular Analysis
    int         moleculestrLength = int(MoleculeName.length());  // String lengeth of molecule
    char        letter;                  // Read MoleculeName from the 1st letter 1by1 -> letter : current letter
    string      nowElement="";           // Store the information of current element
    string      nowNumber="";            // Store the number of element currently reading
    Element     detectedElement;
    
    int _i;
    
    for (int i=0; i<moleculestrLength; i++){
        letter = MoleculeName.at(i);
        if(isupper(letter)){
            // Store nowNu to BcompM when start reading next element
            if (i!=0){
                // cout << nowElement << " " << nowNumber << endl;
                detectedElement = ElementDatabase.correspond(nowElement);
                _i = detectedElement.getNumber();
                if (nowNumber == ""){nowNumber = '1';}
                composition[_i-1] += atoi(nowNumber.c_str());
                
                nowElement.clear();
                nowNumber.clear();
            }
            nowElement.push_back(letter);
        }
        if(islower(letter)){
            nowElement.push_back(letter);
        }
        if(isdigit(letter) || letter=='.'){
            nowNumber.push_back(letter);
        }
        if(letter=='+'){ // do nothing
	}
	if(letter=='('){ break; /*do NOT read the phase information */
	}
    }
    // For the last letter
    detectedElement = ElementDatabase.correspond(nowElement);
    _i = detectedElement.getNumber();
    if (nowNumber == ""){nowNumber = '1';}
    composition[_i-1] = stod(nowNumber.c_str());
    
    /* add up no of elements */
    no_atoms = 0;
    for (int i=0; i<numofElement; i++){ no_atoms += composition[i]; }
    MolecularWeight = 0.0;
    int elist = ElementDatabase.getElementDataSize();
    for (int j=0; j<elist; j++){
	MolecularWeight += composition[j]*ElementDatabase[j].getWeight();
    }
    MolecularWeight *= 1.0e-3;    /* in [kg/mol] */
}
void Molecule::showMoleculeName(){
    cout << "Molecule is " << MoleculeName << endl;
}
Molecule& Molecule::operator=(const Molecule& copy){
    if (this != &copy){
        MoleculeName    = copy.MoleculeName;
        composition     = copy.composition;
        MolecularWeight = copy.MolecularWeight;
    }
    return *this;
}


/*--------------------------------------------------------------------------
 class: ElementData
 ---------------------------------------------------------------------------*/
Element ElementData::operator[](const int i){
    return ElementDatabase[i];
}
Element ElementData::correspond(string name){
    for (int i=0; i<totalNumofElement; i++){
        if(ElementDatabase[i].isthisElement(name)){return ElementDatabase[i];}
    }
    cout << "No corresponding element was found. -> Review your file or Element list : " << name << endl;
    exit(7);
}
ElementData::~ElementData(){
    ElementDatabase.clear();
}

/*--------------------------------------------------------------------------
 class: MoleculeData
 ---------------------------------------------------------------------------*/
MoleculeData::MoleculeData(string filename){
    // Open File
    ifstream fMo(filename);
    if (!fMo){ cout << "Error: Failed to open file (Resulting Molecule Data)" << endl;}
    else{      }//cout << "File Opened (Molecule Data : " << filename << " )" << endl;}
    
    string s_name;
    while(!fMo.eof()){
        fMo >> s_name;
        Molecule data_molecule(s_name);
        MoleculeDatabase.push_back(data_molecule);
    }
    fMo.close();
    
    ElementData ElementDatabase;
    
    numofElement  = ElementDatabase.getElementDataSize();
    numofMolecule = MoleculeDatabase.size();
}


/*--------------------------------------------------------------------------
 // Reference Database
 ---------------------------------------------------------------------------*/
ElementData::ElementData(){
    Element data;
    
    // Hydrogen
    data.setNumber(1);    data.setName("H");
    data.setWeight(1.008);
    ElementDatabase.push_back(data);
    
    // Helium
    data.setNumber(2);    data.setName("He");
    data.setWeight(4.0026);
    ElementDatabase.push_back(data);
    
    // Lithium
    data.setNumber(3);    data.setName("Li");
    data.setWeight(6.941);
    ElementDatabase.push_back(data);
    
    // Beryllium
    data.setNumber(4);    data.setName("Be");
    data.setWeight(9.012);
    ElementDatabase.push_back(data);
    
    // Boron
    data.setNumber(5);    data.setName("B");
    data.setWeight(10.811);
    ElementDatabase.push_back(data);
    
    // Carbon
    data.setNumber(6);    data.setName("C");
    data.setWeight(12.0107);
    ElementDatabase.push_back(data);
    
    // Nitrogen
    data.setNumber(7);    data.setName("N");
    data.setWeight(14.0067);
    ElementDatabase.push_back(data);
    
    // Oxygen
    data.setNumber(8);    data.setName("O");
    data.setWeight(15.9994);
    ElementDatabase.push_back(data);
    
    // Fluorine
    data.setNumber(9);    data.setName("F");
    data.setWeight(18.9984);
    ElementDatabase.push_back(data);
    
    // Neon
    data.setNumber(10);    data.setName("Ne");
    data.setWeight(20.1797);
    ElementDatabase.push_back(data);
    
    // Sodium
    data.setNumber(11);    data.setName("Na");
    data.setWeight(22.9898);
    ElementDatabase.push_back(data);
    
    // Magnesium
    data.setNumber(12);    data.setName("Mg");
    data.setWeight(24.3050);
    ElementDatabase.push_back(data);
    
    // Aluminium
    data.setNumber(13);    data.setName("Al");
    data.setWeight(26.9815);
    ElementDatabase.push_back(data);
    
    // Silicon
    data.setNumber(14);    data.setName("Si");
    data.setWeight(28.0855);
    ElementDatabase.push_back(data);
    
    // Phosphorus
    data.setNumber(15);    data.setName("P");
    data.setWeight(30.9738);
    ElementDatabase.push_back(data);
    
    // Sulfur
    data.setNumber(16);    data.setName("S");
    data.setWeight(32.0650);
    ElementDatabase.push_back(data);
    
    // Chlorine
    data.setNumber(17);    data.setName("Cl");
    data.setWeight(35.453);
    ElementDatabase.push_back(data);
    
    // Argon
    data.setNumber(18);    data.setName("Ar");
    data.setWeight(39.948);
    ElementDatabase.push_back(data);
    
    // Potassium
    data.setNumber(19);    data.setName("K");
    data.setWeight(39.0983);
    ElementDatabase.push_back(data);
    
    // Calcium
    data.setNumber(20);    data.setName("Ca");
    data.setWeight(40.078);
    ElementDatabase.push_back(data);
    
    // Scandium
    data.setNumber(21);    data.setName("Sc");
    data.setWeight(44.9559);
    ElementDatabase.push_back(data);
    
    // Titanium
    data.setNumber(22);    data.setName("Ti");
    data.setWeight(47.867);
    ElementDatabase.push_back(data);
    
    // Vanadium
    data.setNumber(23);    data.setName("V");
    data.setWeight(50.9415);
    ElementDatabase.push_back(data);
    
    // Chromium
    data.setNumber(24);    data.setName("Cr");
    data.setWeight(51.9961);
    ElementDatabase.push_back(data);
    
    // Manganese
    data.setNumber(25);    data.setName("Mn");
    data.setWeight(54.9380);
    ElementDatabase.push_back(data);
    
    // Iron
    data.setNumber(26);    data.setName("Fe");
    data.setWeight(55.845);
    ElementDatabase.push_back(data);
    
    // Cobalt
    data.setNumber(27);    data.setName("Co");
    data.setWeight(58.9332);
    ElementDatabase.push_back(data);
    
    // Nickel
    data.setNumber(28);    data.setName("Ni");
    data.setWeight(58.6934);
    ElementDatabase.push_back(data);
    
    // Copper
    data.setNumber(29);    data.setName("Cu");
    data.setWeight(63.546);
    ElementDatabase.push_back(data);
    
    // Zinc
    data.setNumber(30);    data.setName("Zn");
    data.setWeight(65.38);
    ElementDatabase.push_back(data);
    
    totalNumofElement = ElementDatabase.size();
}
