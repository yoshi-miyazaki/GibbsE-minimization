//  massb_oxide.cpp
//
//  Created by Yoshi Miyazaki on 2017/3/7.
//  Copyright (c) 2017 Yoshi. All rights reserved.

#include "massb.h"

massBalance::massBalance(MoleculeData_G& moleculelist){    
    create_massBalance(moleculelist);
}
void massBalance::create_massBalance(MoleculeData_G& list){

    int m = list.getnumofMolecule();
    tensor1d<int> check(0,numofElement);
    /* show whether i'th element is used in the system or not */
    ElementData elementdatabase;
    for (int i=0; i<numofElement; i++){
        for (int j=0; j<m; j++){
            if(list[j].getComposition(i) != 0){check[i] = 1;}
        }
    }
    
    /* d: num of elements used in this system = sum of vec check */
    int d = check.sum();
    massm.resize(0,m,d);  atomic.resize(0,d);
    
    /* Store molecule information into a matrix */
    int k = 0;  /* column # in the resized matrix */
    for (int i=0; i<numofElement ; i++){
        if (check[i] == 1){ /* removed elements that are not used in any species */
            for (int j=0; j<m; j++){ massm[j][k] = list[j].getComposition(i); }
            atomic[k] = i+1;  /* store the actual atmoic # (not 0 but 1 for H) */
            k++;
        }
    }
    
    /*cout << "Mass balance (massb_oxide.cpp)" << endl;	*/
}
void massBalance::create_nullspace(MoleculeData_G& moleculelist){
    m = moleculelist.getnumofMolecule();
    d = m - getdim(moleculelist);        // dim of null-space
    nspace.resize(0.0,m,d);
    which.resize(m);

    int i_dim = 0;   /* # of nullspace vector created */
    for (int i=0; i<m; i++){
        Molecule_Gibbs mol = moleculelist[i];
        int i_simple = issimple(mol);
        if (i_simple == 0){ /* if compound */

	    which[i] = i_dim;  

	    /* Create null space vector
               componud +1 and balance w/ simple
               ex. Al2O3 +1
                   Al    -2
	           O2    -1.5 */
            nspace[i][i_dim] += 1;
            
            for (int j=0; j<numofElement; j++){
                if (mol.getComposition(j) != 0){
                    for (int k=0; k<m; k++){
                        if(moleculelist[k].getMoleculeName() == simple(j)){
                            nspace[k][i_dim] -= mol.getComposition(j)/simple_index(j);
                        }
                    }
                }
            }
                        
            // finish creating 1 nullspace vector -> increase i_dim
            i_dim++;
            if (i_dim > d){cout << "Molecule list incomplete. Include corresponding simple substance for all the element."  << " " << i_dim << "  , d: " << d << " , m: " << m << endl; exit(25);}
        }
	else{
	    which[i] = -1;
	}
    }

    cout << "Null space" << endl;
    for (int i=0; i<m; i++){
        for (int j=0; j<d; j++){
            cout << setw(6) << nspace[i][j] << " " ;
        }
        cout << endl;
    }
}

/*---------------------------------------------------
  Count the number of existing mols
  ex.) If H2, O2, H2O, CO are in ths system
       ---> H, C, O so return "3"
 ---------------------------------------------------*/
int massBalance::getdim(MoleculeData_G& moleculelist){
    
    tensor1d<double> nonzero(0.0,numofElement);    /* i'th element exist in this system */
                                                   /* ---> nonzero != 0                 */

    int m = moleculelist.getnumofMolecule();
    for (int i=0; i<m; i++){
        tensor1d<double> temp(numofElement);
        for (int j=0; j<numofElement; j++){
            temp[j] = moleculelist[i].getComposition(j);
        }
        nonzero += temp;
    }
    
    int count = 0;
    for (int j=0; j<numofElement; j++){
        if (nonzero[j] != 0){ count++; }
    }
    return count;
}

int massBalance::issimple(Molecule& mol){
    string _s = mol.getMoleculeName();
    
    if (_s == "H2"){return 1;}
    if (_s == "He"){return 1;}
    if (_s == "Li"){return 1;}
    if (_s == "Be"){return 1;}
    if (_s == "B") {return 1;}
    if (_s == "C") {return 1;}
    if (_s == "N2"){return 1;}
    if (_s == "O2"){return 1;}
    if (_s == "F2"){return 1;}
    if (_s == "Ne"){return 1;}
    if (_s == "Na"){return 1;}
    if (_s == "Mg"){return 1;}
    if (_s == "Al"){return 1;}
    if (_s == "Si"){return 1;}
    if (_s == "P") {return 1;}
    if (_s == "S") {return 1;}
    if (_s == "Cl2"){return 1;}
    if (_s == "Ar"){return 1;}
    if (_s == "K") {return 1;}
    if (_s == "Ca"){return 1;}
    if (_s == "Sc"){return 1;}
    if (_s == "Ti"){return 1;}
    if (_s == "V") {return 1;}
    if (_s == "Cr"){return 1;}
    if (_s == "Mn"){return 1;}
    if (_s == "Fe(g)"){return 1;}
    if (_s == "Co"){return 1;}
    if (_s == "Ni"){return 1;}
    if (_s == "Cu"){return 1;}
    if (_s == "Zn"){return 1;}
    
    // else
    return 0;
}

int massBalance::simple_index(const int k){
    if (k == 0){return 2;}  // H2
    if (k == 6){return 2;}  // N2
    if (k == 7){return 2;}  // O2
    if (k == 8){return 2;}  // F2
    if (k ==16){return 2;}  // Cl2
    
    return 1;
}

// Return k'th element's simple substance
string massBalance::simple(const int k){
    if (k == 0){return "H2";}
    if (k == 1){return "He";}
    if (k == 2){return "Li";}
    if (k == 3){return "Be";}
    if (k == 4){return "B";}
    if (k == 5){return "C";}
    if (k == 6){return "N2";}
    if (k == 7){return "O2";}
    if (k == 8){return "F2";}
    if (k == 9){return "Ne";}
    if (k ==10){return "Na";}
    if (k ==11){return "Mg";}
    if (k ==12){return "Al";}
    if (k ==13){return "Si";}
    if (k ==14){return "P";}
    if (k ==15){return "S";}
    if (k ==16){return "Cl2";}
    if (k ==17){return "Ar";}
    if (k ==18){return "K";}
    if (k ==19){return "Ca";}
    if (k ==20){return "Sc";}
    if (k ==21){return "Ti";}
    if (k ==22){return "V";}
    if (k ==23){return "Cr";}
    if (k ==24){return "Mn";}
    if (k ==25){return "Fe";}
    if (k ==26){return "Co";}
    if (k ==27){return "Ni";}
    if (k ==28){return "Cu";}
    if (k ==29){return "Zn";}
    
    return "H2";
}
