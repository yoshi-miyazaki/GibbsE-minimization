/*
  initcomp.cpp
  Gibbs
  
  Created by Yoshi Miyazaki on 2015/04/11.
  Revised on Feb 25, 2016
  Copyright (c) 2015 Yoshi. All rights reserved.
*/

#include "initcomp.h"

initcomp::initcomp(string init_comp, MoleculeData_G& moleculelist){
    // cout << endl << "--- [ initial composition for minimization ] ---" << endl;
    list_species = moleculelist.getlist();   /* list of molecules stored in initcomp */
    
    // Open input file (init_comp)
    // and save information in 2 different objects
    //       1. MoleculeData (molecule information)
    //       2. initial_mol  (mol # of each molecule)
    // done by function read_initcomp

    /* list of molecule included in Composition.txt */
    MoleculeData     initmoleculelist;
    initmoleculelist = read_initcomp(init_comp);
    
    /* match "initmoleculelist" and "moleculelist" */
    int numofconsideringmol = moleculelist.getnumofMolecule();
    output_mol.resize(0.0,numofconsideringmol);
    
    // cout << "initial composition file: " << init_comp << endl;
    for (int i=0; i<numofcomposition; i++){
        Molecule _m = initmoleculelist[i];
        // cout << _m.getMoleculeName() << " - " << i << " / " << (numofcomposition-1) << endl;
        for (int j=0; j<numofconsideringmol; j++){
            if (_m == moleculelist[j]){
                string _s = _m.getMoleculeName();
                
                // if molecules are...
                if (_s == "H"){ for (int k=0; k<numofconsideringmol; k++){
                        if (moleculelist[k].getMoleculeName() == "H2") {
                            output_mol[k] += initial_mol[i]/2;}}}
                else if (_s == "N"){ for (int k=0; k<numofconsideringmol; k++){
                        if (moleculelist[k].getMoleculeName() == "N2") {
                            output_mol[k] += initial_mol[i]/2;}}}
                else if (_s == "O"){ for (int k=0; k<numofconsideringmol; k++){
                        if (moleculelist[k].getMoleculeName() == "O2") {
                            output_mol[k] += initial_mol[i]/2;}}}
                else if (_s == "F"){ for (int k=0; k<numofconsideringmol; k++){
                        if (moleculelist[k].getMoleculeName() == "F2") {
                            output_mol[k] += initial_mol[i]/2;}}}
                else if (_s == "Cl"){ for (int k=0; k<numofconsideringmol; k++){
                        if (moleculelist[k].getMoleculeName() == "Cl2") {
                            output_mol[k] += initial_mol[i]/2;}}}
                else{   output_mol[j] += initial_mol[i];}
                break;
            }
        }
    }
    
    /* output to console and Check 
    for (int j=0; j<numofconsideringmol; j++){
        cout << setw(10) << moleculelist[j].getMoleculeName() << " = " << output_mol[j] << endl;
	}*/
        
    /* element composition ... add up each molecule information */
    element_mol.resize(0.0,numofElement);
    for (int p=0; p<numofElement; p++){
        element_mol[p] = 0;
        for (int i=0; i<numofcomposition; i++){
            element_mol[p] += initial_mol[i] * initmoleculelist[i].getComposition(p);
        }
    }
}

// Separate the file into moleculename-file and num_of_mol data
MoleculeData initcomp::read_initcomp(string init_comp){
    // open File
    ifstream fcomp(init_comp);
    if(!fcomp){ cout << "Error: Failed to open initcomp file: " << init_comp << endl;}
    else{   }//     cout << "File Opened (Initial Composition)" << endl;}
    
    // store string(moleculename) to s_name
    // & double(mol#) to n_comp
    string   s_name, p_name;
    double   n_comp;
    string   outputfile="moleculelist.txt";
    ofstream fout(outputfile);
    
    // cout << "reading. " << endl;
    numofcomposition = 0; /* initialize */
    while(!fcomp.eof()){
        fcomp >> s_name >> n_comp;  /* read each line */
        if     (s_name == "" or s_name == p_name){ cout << "initcom.cpp -- read twice..." << s_name << endl; break; }
        else{ numofcomposition++; }
        
        // store molecule name list
        // do not store mol# in the file.
        fout  << endl << s_name;            // Read the molecule name and store it in a new file
        // store "<< endl" first so that "\n" won't be inserted in the end of the file
        //              --> that would cause to read the last line of the data twice.
        
        // store mol #
        Vector1d<double> temp(numofcomposition);
        if (numofcomposition != 1){         // Copy the previous data
            for (int i=0; i<(numofcomposition-1); i++){ temp[i] = initial_mol[i]; }
        }
        temp[numofcomposition-1] = n_comp;  // add the data of newest molecule
        initial_mol = temp;
        
        /* avoid reading the last line twice */
        p_name = s_name;
    }
    fout.close();
    
    // prepare file ONLY w/ moleculename and create a object MoleculeData w/ that file
    MoleculeData initmoleculelist(outputfile);
    return initmoleculelist;
}
