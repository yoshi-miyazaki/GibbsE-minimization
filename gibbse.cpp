//
//  gibbse.cpp
//
//  Created by Yoshi Miyazaki on 2015/04/10.
//  Copyright (c) 2015 Yoshi. All rights reserved.
//

#include "gibbse.h"

extern const int numofElement;


/*--------------------------------------------------------------------------
  class: MoleculeData_G
 ---------------------------------------------------------------------------*/

MoleculeData_G::MoleculeData_G(string filename, double T, double P){
    // Open File
    clock_t  begin = clock();
    ifstream fMo(filename);
    if (!fMo){ cout << "Error in MoleculeData_G: Failed to open file (Resulting Molecule Data)" << endl;}
    else{ } //cout << "File Opened (Molecule Data : " << filename << " )" << endl;}
    
    // Read line one by one
    while(!fMo.eof()){
        string s_name;
        fMo >> s_name;
        if (s_name == ""){break;}
        Molecule_Gibbs data_molecule(s_name, T, P);
        database.push_back(data_molecule);
    }
    fMo.close();
    
    numofMolecule = database.size();     // # of molecule
    list_species  = filename;

    molarmass.resize(numofMolecule);
    for (int i=0; i<(int)numofMolecule; i++){ molarmass[i] = database[i].getWeight(); }

    clock_t end = clock();
    // cout << "time to create GibbseE database: " << double(end-begin)/CLOCKS_PER_SEC << " - " << numofMolecule << " @T=" << T << endl;
}
MoleculeData_G& MoleculeData_G::operator=(const MoleculeData_G& copy){
    if (this != &copy){
        numofMolecule = copy.numofMolecule;
        database      = copy.database;
    }
    return *this;
}
int  MoleculeData_G::intspec(string _s){ /* return index of species _s */
    int index = -1;
    for (int i=0; i<(int)numofMolecule; i++){
        if (_s == database[i].getMoleculeName()){ index = i; break; }
    }
    if (index == -1){ } // cout << "WARNING: no species " << _s << endl; }
    return index;
}
void MoleculeData_G::erase(int i){
    database.erase(database.begin()+i);
    numofMolecule = database.size();
    cout << "deleted no. " << i << " , now list-size: " << numofMolecule << endl;
}
void MoleculeData_G::rm_nonexisting(Vector1d<double>& ninit){
    /* remove species that won't appear in the given set of elemental composition
       e.g.) FeO should be removed if there are no Fe in the system... 
    
       so, check for each considering species (int j, j'th molecule[p] != 0)
       whether it includes any non-existing elements (int p, element_mole[p] == 0) */

    /* element composition ... add up each molecule information */
    Vector1d<double> element_mole(0.0,numofElement);
    for (int p=0; p<numofElement; p++){
        for (int i=0; i<(int)numofMolecule; i++){
            element_mole[p] += ninit[i] * database[i].getComposition(p);
        }
    }
    
    /* decide on the species that should be removed. */
    Vector1d<int> rm_species(0,numofElement);
    for (int i=0; i<(int)numofMolecule; i++){
        for (int p=0; p<numofElement; p++){
            int element_p = database[i].getComposition(p);
            if (element_mole[p] == 0 && element_p != 0){ rm_species[i] = 1; break; }
        }
    }
    
    /* adjust starting vector ninit as well */
    for (int j=(numofMolecule-1); j>=0; j--){ /* remove from species w/ larger index */
        if (rm_species[j] == 1){
            numofMolecule--;                    /* decrease list size */
            database.erase(database.begin()+j); /* delete Molecule    */
            ninit.erase(j); }                   /* adjust n_init      */
    }
    // cout << "no. of Molecule = " << numofMolecule << endl;
    // for (int j=0; j<numofMolecule; j++){ cout << database[j].getMoleculeName() << endl; }
}
void MoleculeData_G::rm_nonexisting(tensor1d<double>& ninit){
    /* remove species that won't appear in the given set of elemental composition
       e.g.) FeO should be removed if there are no Fe in the system... 
    
       so, check for each considering species (int j, j'th molecule[p] != 0)
       whether it includes any non-existing elements (int p, element_mole[p] == 0) */

    /* element composition ... add up each molecule information */
    Vector1d<double> element_mole(0.0,numofElement);
    for (int p=0; p<numofElement; p++){
        for (int i=0; i<(int)numofMolecule; i++){
            element_mole[p] += ninit[i] * database[i].getComposition(p);
        }
    }
    
    
    /* decide on the species that should be removed. */
    Vector1d<int> rm_species(0,numofMolecule);
    for (int i=0; i<(int)numofMolecule; i++){
        for (int p=0; p<(int)numofElement; p++){
            int element_p = database[i].getComposition(p);
            if (element_mole[p] == 0 && element_p != 0){ rm_species[i] = 1;  break; }
        }
    }
    
    /* adjust starting vector ninit as well */
    for (int j=(numofMolecule-1); j>=0; j--){ /* remove from species w/ larger index */
        if (rm_species[j] == 1){
            numofMolecule--;                    /* decrease list size */
            database.erase(database.begin()+j); /* delete Molecule    */
            ninit.erase(j); }                   /* adjust n_init      */
    }
    // cout << "no. of Molecule = " << numofMolecule << endl;
    // for (int j=0; j<numofMolecule; j++){ cout << database[j].getMoleculeName() << endl; }
}
/*--------------------------------------------------------------------------
  class: Molecule_Gibbs
 ---------------------------------------------------------------------------*/
Molecule_Gibbs::Molecule_Gibbs(string _s, double Ti, double Pi) : Molecule(_s){    
    // Search for Gibbs Energy
    T = Ti;  P = Pi;                 /* store T&P in this object  */
    gibbsE = data_search(_s);  /* gibbsE is in [J] NOT [kJ] */
    
    // convert Gibbs energy into given condition
    double Pbar = P/1e5;   /* convert to bar */
    if (phase == 0){gibbsE /= (R*T); gibbsE += log(Pbar);} // gas    phase
    if (phase == 1){gibbsE /= (R*T);                     } // liquid phase
    if (phase == 2){gibbsE /= (R*T);                     } // solid  phase
}
double Molecule_Gibbs::data_search(string s_mol){
    /* read Gibbs data from various data format
     ... 1. first from Stixrude */
    string file1  = "/Users/Yoshi/Documents/OneDrive/OneDrive/Programs/C++/GibbsFE_minimization/Gibbs_data/Gibbs_Stixrude.txt";
    // string file1 = "/nfs4mounts/home/ym267/GibbsFE_minimization/Gibbs_data/Gibbs_Stixrude.txt";
    
    /* first check Stixrude database */
    string s_name;   char c_phase;
    double Gf = 0.0;
    ifstream fStixrude(file1);
    if (!fStixrude){ cout << "gibbse.cpp- No Stixrude database: " << file1 << endl; exit(7);}
    else{
        double F0, V0, debye0, gma0, q0;
        double dG;
        int    read = 0;
        while(!fStixrude.eof()){
            /* Stixrude database is written in the order of:
               moleculename(phase)  s/l/v-phase  Helmholtz[kJ/mol]  molar-volume[cm3/mol]
               bulk-modulus[GPa]  K0'  Debye-T[K]  gamma-0  q-0*/
            fStixrude >> s_name >> c_phase >> F0 >> V0 >> K0 >> K0d >> debye0 >> gma0 >> q0 >> dG;
            // cout << "reading... " << s_name << endl;
	    
            if (s_name == s_mol){read = 1; break;}
        }
        if (read == 1){
            /* if found on the list, put into equation of state */
            V0 *= 1.0e-6;  F0 *= 1.0e3;  K0 *= 1.0e9; /* convert to SI unit */
            Gf = GibbsE_Stixrude(F0, V0, K0, K0d, debye0, gma0, q0) +dG;
            set_phase(c_phase);
            fStixrude.close();
            return Gf;
        }
    }
    
    /* next check JANAF database
       ... open database (in ./Gibbs_data file. extension:.txt) */
    string fileJ = "/Users/Yoshi/Documents/OneDrive/OneDrive/Programs/C++/GibbsFE_minimization/Gibbs_data/" + s_mol + ".txt";
    // string fileJ = "/nfs4mounts/home/ym267/GibbsFE_minimization/Gibbs_data/" + s_mol + ".txt";
    
    ifstream fGibbs(fileJ);
    /* check the existence of file. If no file found, move on. */
    if (fGibbs){
        tie(Gf, c_phase) = GibbsE_JANAF(fGibbs, s_mol);
        set_phase(c_phase);
        return Gf;
    }
    
    /* for SiO2(l) use the original formulation */
    if (s_mol == "SiO2(l)"){
        double F0 = 0.0, V0 = 24.5e-6, S0 = 41.463;
        // double K0 = 16.24e9, K0d = 5.0, debye0 = 400.0, gma0 = 0.78, q0 = 1.0;
        double K0 = 25.1e9, K0d = 4.41, alpha0 = 5.82e-6, C0 = 2.13, C1 = 0.0553, d0 = 1;
        Gf = GibbsE_EoS(F0, V0, S0, K0, K0d, alpha0, C0, C1, d0);
	
        double Ttmp = T, Ptmp = P;  T = 1943.; P = 1.e5;
        Molecule_Gibbs qtz("SiO2(q)",T,P);
        double Gqtz = qtz.getGibbsE()*(R*T);
        double Gliq = GibbsE_EoS(F0, V0, S0, K0, K0d, alpha0, C0, C1, d0);
        T = Ttmp, P = Ptmp;
        Gf -= (Gliq-Gqtz);
	
        set_phase('l');
        return Gf;
    }
    
    /* if not found, look for liquid-database file */
    string file2 = "/Users/Yoshi/Documents/OneDrive/OneDrive/Programs/C++/GibbsFE_minimization/Gibbs_data/melting-curve/" + s_mol + ".txt";
    // string file2 = "/nfs4mounts/home/ym267/GibbsFE_minimization/Gibbs_data/melting-curve/" + s_mol + ".txt";
    
    ifstream fmelt(file2);
    if (fmelt){
        Gf = GibbsE_melting(fmelt);
        c_phase = 'l';  set_phase(c_phase);
        fmelt.close();
        return Gf;
    }
    
    cout << "(No GibbsFE database for molecule " << s_mol << ".)" << endl;  exit(9);
    return Gf;
}

/* ---------- basic functions. used in all 3 ----------*/
void Molecule_Gibbs::set_phase(char c_phase){
    if(c_phase=='s'){phase = 2;}
    if(c_phase=='l'){phase = 1;}
    if(c_phase=='g'){phase = 0;}    
}
double Molecule_Gibbs::internal_division(double T, double T1, double T2, double d1, double d2){
    /*--------------------------------------------------------------
    // Internal & External Division
    //                                              x (d2)
    //
    //                         x (d1)
    //            (d)o
    //         <- calculate this value
    //       --------|---------|--------------------|
    //               T        T1                   T2
    --------------------------------------------------------------*/
    double d = d1 + (d2-d1)*(T-T1)/(T2-T1);
    return d;
}

/* ---------- reading from Stixrude's database ----------*/
double Molecule_Gibbs::GibbsE_Stixrude(double F0, double V0, double K0, double K0d, double debye0, double gma0, double q0){
    /* calculate Gibbs free energy using EoS by Stixrude (2011) */
    double rho0 = 1/V0;
    double rhop = rho0*pow(1.0+K0d*P/K0, 1/K0d);  /* molar density at P=P [m3/mol] */
    double fi   = 0.5*(pow(rhop/rho0,2.0/3) -1);  /* strain */
    double F = F0 + V0*(4.5*K0*(fi*fi)*(1+(K0d-4)*fi));
    
    double dDebye = T*debyeF_func(T,fi,debye0,gma0,q0) - Tref*debyeF_func(Tref,fi,debye0,gma0,q0);
    int    nmol = no_atoms;
    double FdT = 9*nmol*kB*NA*dDebye;
    double Gf = F + FdT + P/rhop;  /* Helmholtz cold, thermal, and P*V for Gibbs conversion */
    //cout << MoleculeName << " - Gf: " << Gf << " - fi:" << fi << " - rhop: " << rhop << " do: " << pow(1.0+K0d*P/K0, 1/K0d) << ", P: " << P <<endl;
    
    return Gf;
}
double Molecule_Gibbs::GibbsE_EoS(double F0, double V0, double S0, double K0, double K0d, double alpha0, double C0, double C1, double d0){
    /* calculate Gibbs free energy using EoS by Stixrude (2011) */
    double VP = V0*K0/(K0d-1)*(pow(1.0+K0d*P/K0, 1-1/K0d)-1);
    double ST = S0*(T-Tref) +C0*(T*log(T/Tref) +(T-Tref)) +0.5*C1*square(T-Tref);
    double S_P = int_dSdP(V0,S0,K0,K0d,alpha0,d0);
    
    double alpha = alpha0*exp(-1.*d0/K0*P);
    double Gf = F0 + VP - ST - S_P/alpha*(exp(alpha*(T-Tref))-1);
    /* Helmholtz + cold compression + thermal effect + P correction to thermal effect */
    //cout << MoleculeName << " - Gf: " << Gf << " - fi:" << fi << " - rhop: " << rhop << " do: " << pow(1.0+K0d*P/K0, 1/K0d) << ", P: " << P <<endl;
    // cout << "S0 ~ " << S0 << " , correction: " << S_P << " V0: " << V0 << " , K0: " << K0 << " , T0 = " << T0 << " , T = " << T << " , d0 ~ " << d0 << endl;
    
    return Gf;
}
double Molecule_Gibbs::int_dSdP(double V0, double S0,  double K0, double K0d, double alpha0, double d0){
    /* integrate dS/dP from Patm to P */
    int N = 100;    double dP = (P-Patm)/double(N-1);
    Vector1d<double>  Pd(0.0,N);
    
    double S = 0;
    for (int i=0; i<N; i++){
	Pd[i] = Patm + dP*(double)i;
        S -= alpha0*V0*exp(-1.*d0/K0*Pd[i])/pow(1+K0d*Pd[i]/K0,1/K0d)*dP;
	// cout << " - P = " << Pd[i]/1e9 << " al/al0 = " << exp(-1.*d0/K0*Pd[i]) << " - " << pow(1+K0d*Pd[i]/K0,1/K0d) << " - " << dP << " total = " << exp(-1.*d0/K0*Pd[i])/pow(1+K0d*Pd[i]/K0,1/K0d)*dP << endl;
    }
    return S;
}
double Molecule_Gibbs::GibbsE_melting(ifstream& fmelt){
    /* read melting temperature from fmelt, and convert into Gibbs free energy data */
    double Gf;
    string sline; /* store each line */

    double V0, DS0, A;   /* read the 1st line, molar volume, entropy change at 0Pa,
			             d(dS0)/dP are stored */
    fmelt >> V0 >> DS0 >> A;
    V0 *= 1.0e-6;    /* cm3/mol to m3/mol */
    
    double Tc, Tb = 0, Ta = 0;
    double Pc, Pb = 0, Pa = 0;
    string s_prev, s_solid;
    while(getline(fmelt, sline)){
	/* in each (l) file, data of melting curve is stored
	   ... pressure [GPa], temperature [K], and the corresponding solid phase */
	if (sline.empty()){ continue; }
	istringstream smelt(sline);
	smelt >> Pc >> Tc >> s_solid;
	Pc *= 1.0e9; /* convert from GPa to Pa */
	
	if (P < Pc){
	    // cout << "D - Pc: " << Pc/1e9 << " GPa , " << Pb/1e9 << endl;
	    if (Pb == 0.0){
		/* if the 1st pressure in the list is lower than the given pressure,
		 extrapolate using the 1st and the 2nd line in the melting curve data */
		Tb = Tc;  Pb = Pc;  s_prev = s_solid;
		getline(fmelt, sline);
		istringstream smelt(sline);
		smelt >> Pc >> Tc >> s_solid;
		Pc *= 1.0e9; /* convert from GPa to Pa */
		// cout << "Pc: " << Pc/1e9 << " GPa , " << Pb/1e9 << endl;
	    }
	    break;
	}
	/* for last row, in case of double reading */
	if (Pc == Pb && Pc != 0){ Tb = Ta;  Pb = Pa; break; } 
	
	/* prepare to read next row */
	Ta = Tb;  Pa = Pb;
	Tb = Tc;  Pb = Pc;    s_prev = s_solid;
    }
    if (P > Pc){ Tb = Ta; Pb = Pa; } /* when T is over the data range,
					extrapolate using last two lines  */
    
    /* calculate melting tempeature at T=T from (Tprev, Pprev) and (Tm,Pm) */
    double Tmelt = internal_division(P, Pb, Pc, Tb, Tc);
    // cout << "* - calc melting of " <<  s_prev << endl;
    // cout << "Melting T at: " << setprecision(5) << Tmelt << " K at P = " << P << "[K]. b/w " << Tb << " and " << Tc << " K, P =" << Pb/1e9 << " and " << Pc/1e9 << "GPa. " << endl;
    
    /* Gibbs free energy should be the same b/w solid and liquid at melting temperature.
       Gibbs energy of liquid would be calculated by adding 
       the difference of Gibbs energy from the melting temperature */
    Molecule_Gibbs mole_solid(s_prev, T, P);
    
    /* calc liquid entropy as Sliq = Ssol + dS
       ... use the fact that dS does not change much w/ temperature.
           we assume dS as a simple function of pressure
	   dS = dS0 [1 - A*(1 - exp(-P/P0))] 
    
       ... in dataset, dS0 and A must be specified. 
           P0 is set to 100 GPa in our model. */
    double P0 = 50e9;
    double DS = DS0*(1 - A*(1 - exp(-1.*P/P0)));
    
    double Gsol_Tmelt = GibbsE_solid(s_prev,Tmelt,P);
    double Gsol_T     = GibbsE_solid(s_prev,T,    P);	
    
    double negSdT_sol = Gsol_T - Gsol_Tmelt;
    double negSdT_liq = negSdT_sol - DS*(T-Tmelt);
    
    /* Use the fact that G(sol)(Tmelt,P) = G(liq)(Tmelt,P)
       
       Tmelt |----------#    *
             |       *  |
             |   *      |
           T |*---------#  ... adjust T along pressure P
             |          |       => adjust T using Sliq
             |          |
             |__________|_______
                    Pb  P   Pc
    */
	
    Gf = Gsol_Tmelt + negSdT_liq;
    // cout << "   G of : sol: " << Gsol_T/(R*T) << " - liq: " << Gf/(R*T) << endl << endl;
    
    return Gf;
}
double Molecule_Gibbs::GibbsE_solid(string s_prev, double Tmelt, double Pmelt){
    /* calculate GFE of s_prev (solid) at T=Tmelt, P=Pmelt */
    Molecule_Gibbs mole_solid(s_prev, Tmelt, Pmelt);
    double Gsolid = mole_solid.getGibbsE();
    
    Gsolid *= (R*Tmelt);  /* convert into [J], the stored value is divided by R*T */
    return Gsolid;
}
double Molecule_Gibbs::int_VdP(double Vl, double K0, double K0d, double P1, double P2){
    /* calculate dG = d(- K/(K0'-1)*V )
       ... this equation is derived from thermodynamic relation: dG/dP = V = V0 (1 + K0d*P/K0)^(-1/K0d) */
    /* V0 here is a few % larger than solid at T=T, P=P1 */
    double V = Vl*pow(1+K0d*P1/K0, 1/K0d); /* convert to P=0, and re-convert to P1 and P2 */

    double V1 = V/pow(1+K0d*P1/K0, 1/K0d), V2 = V/pow(1+K0d*P2/K0, 1/K0d);
    double G1 = V1/(K0d-1)*(K0+K0d*P1),    G2 = V2/(K0d-1)*(K0+K0d*P2);
    
    cout << "volume : " << V1*1e6 << " - dG = " << G2-G1  << " - P1: " << P1/1e9 << " - P2: " << P2/1e9 << " - K0: "  << K0/1e9 << endl;
    cout << "G1: " << G1 << ", G2: " << G2 << " V1: " << V1*1e6 << " , V2: " << V2*1e6 << ", original V: " << V*1e6 << endl;
    return G2-G1;
}
double Molecule_Gibbs::dGdT(string sr, double Tr, double Pr){
    double dT = 1.0e-3;
    Molecule_Gibbs m_base(sr, Tr,    Pr);
    Molecule_Gibbs m_diff(sr, Tr+dT, Pr);
     
    double G_base = m_base.getGibbsE() *(R*Tr);
    double G_diff = m_diff.getGibbsE() *(R*(Tr+dT));
    double dG = G_diff - G_base;
    
    return -1*dG/dT;
}
double Molecule_Gibbs::dGdP(string sr, double Tr, double Pr){
    double dP = 1.0e2;
    Molecule_Gibbs m_base(sr, Tr, Pr);
    Molecule_Gibbs m_diff(sr, Tr, Pr+dP);
     
    double G_base = m_base.getGibbsE() *(R*Tr);
    double G_diff = m_diff.getGibbsE() *(R*Tr);
    double dG = G_diff - G_base;
    
    return dG/dP;
}
double Molecule_Gibbs::debyeF_func(double T, double f, double debye0, double gamma0, double q0){
    /* Use debye temperature to calculate the temperature dependence of Helmholtz free energy */
    /* a1 =  6*gm0, a2 = -6*gm0 +18*gm0^2 -9*q0*gm0
       debye = deb0*sqrt(1 +a1*f +a2*f^2) */
    /* debye function D(x):
       D(x) = 1/x^3 * int[ln(1-exp(-t))*t^2] from 0 to x
       x = debye/T; */
    
    double debye = debye0*sqrt(1 +6*gamma0*f + gamma0*(-6 +18*gamma0 -9*q0)*(f*f));
    double x = debye/T;
    double N = 1000, dt = x/N;
    
    double D = 0.0;
    for(int i=0; i<N; i++){
	double tl = (i)*dt;
	double tr = (i+1)*dt;
	if (tl==0){ tl = tr; }
	
	double fl = log(1-exp(-1.*tl))*(tl*tl);
	double fr = log(1-exp(-1.*tr))*(tr*tr);
	D += 0.5*(fl+fr)*dt;
    }
    
    D /= (x*x*x);
    return D;
}

/* ---------- for JANAF database ---------- */
tuple<double, char> Molecule_Gibbs::GibbsE_JANAF(ifstream& fGibbs, string s_molName){
    /* 1st line of the data file: molecule Name
       Check whether it matches the filename. If did not match, abort execution */
    string    s_readName;
    fGibbs >> s_readName;
    if (s_molName != s_readName){
        cout << "File title and the first line does NOT match. Revise the data." << endl;
        exit(8);
    }
    
    /* 2nd line (T	Cp	S	-[G-H(Tr)]/T	H-H(Tr)     ∆f H	∆f G	log Kf) */
    string _s;
    fGibbs >> _s >> _s >> _s >> _s >> _s >> _s >> _s >> _s;   /* just read */
    
    /* 3rd - last line ... read the data 
       variables to store the data for current line (d_), previous line (l_), and 2nd previous line (l2)
       l_ denotes left */
    double l2temp = 0,  l2Gibbs,  l2aGibbs;  //                    1500(l2temp)     Cp  S  GHr/T  HHr  ∆fH	∆fG(l2Gibbs)	 log Kf
    double l_temp = 0,  l_Gibbs,  l_aGibbs;  //                    1600(l_temp)     Cp  S  GHr/T  HHr  ∆fH	∆fG(l_Gibbs)	 log Kf
    double d_temp = 0,  d_Gibbs,  d_aGibbs;  // Current line -->   1700(d_temp)     Cp  S  GHr/T  HHr  ∆fH	∆fG(d_Gibbs)	 log Kf
    
    double _d, INITIAL = 1e8, formationG = INITIAL;   char c_phase = 'g';
    while (!fGibbs.eof()){
        /* Store the previous data to l2temp, l2Gibbs, l_temp, l_Gibbs
	   1st time : l2temp and l_temp is both 0
	   2nd time : l2temp is 0                                        */
        l2temp  = l_temp;   l_temp  = d_temp;
        l2Gibbs = l_Gibbs;  l_Gibbs = d_Gibbs;
        
        // Read line
        fGibbs >> d_temp >> _d >> _d >> _d >> _d >> _d >> d_Gibbs >> _d >> c_phase;
        if (T < d_temp){
            if (l_temp == 0){  /* CASE 1 */
                /* When the target T is out of data range
		   =  when the first temperature in the data is larger than target temperature
		   -> pursue external division using the second temperature in the data
		   ex.) When the data started from 280K and
		   T     Cp     S       GHr/T   HHr     ∆fH	  ∆fG	  log Kf
		   280   75.56  65.215  70.102  -1.368  -286.41  -240.12  44.796   l
		   the target T was below 280K. */
                
                /* Read the next line and pursue external division */
                l_temp  = d_temp;
                l_Gibbs = d_Gibbs;
                fGibbs >> d_temp >> _d >> _d >> _d >> _d >> _d >> d_Gibbs >> _d >> c_phase;
                // cerr << "WARNING: " << T << " is out of thermodynamics table data range (molecule " << s_molName << ".)" << endl;
		
                /* Internal and external division share the same equation */
                formationG = internal_division(T,l_temp,d_temp,l_Gibbs,d_Gibbs);
                break;
            }
            else {  /* CASE 2 */
                /* Internal division
		   ex.) When target T : 1650K  was covered by 1600K and 1700K data. */
                formationG = internal_division(T,l_temp,d_temp,l_Gibbs,d_Gibbs);
                break;
            }
        }
        /* CASE 3 : When the data exists for target T. */
        else if (d_temp == T){
            formationG = d_Gibbs;
            break;
        }
    }
    
    /* CASE 4 : When the target T is higher than the data range
	 Use 2 last data to extrapolate
	 external division              */
    if (formationG == INITIAL){
	/* note that it Does NOT work when formation Gibbs E = INITIAL, which would NOT happen */
        // cout << "WARNING Upper : temperature is out of thermodynamics table data range (molecule " << s_molName << ".)" <<  endl;
	
        // Use l2temp and d_temp to extrapolate (-> l_temp is overwritten by d_temp at this point)
        formationG = internal_division(T,l2temp,d_temp,l2Gibbs,d_Gibbs);
    }
    
    // Set the phase
    if (c_phase == 'g'){phase = 0;}
    if (c_phase == 'l'){phase = 1;}
    if (c_phase == 's'){phase = 2;}
    
    fGibbs.close();
    return forward_as_tuple(formationG*1e3, c_phase);
}


/* ---------- unused functions ---------- */
double Molecule_Gibbs::dH1(double T1, double T2, double A, double B, double C, double D){
    double dH;
    B *= 1e-3; C *= 1e5; D *= 1e-6;
    dH = A*(T2-T1) + B/2*(T2*T2-T1*T1) - C*(1/T2-1/T1) + D/3*(T2*T2*T2-T1*T1*T1);
    return dH;
}
double Molecule_Gibbs::dS1(double T1, double T2, double A, double B, double C, double D){
    double dS;
    B *= 1e-3; C *= 1e5; D *= 1e-6;
    dS = A*(log(T2)-log(T1)) + B*(T2-T1) - C/2*(1/(T2*T2)-1/(T1*T1)) + D/2*(T2*T2-T1*T1);
    return dS;
}
