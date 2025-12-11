//
//  gibbse.cpp
//
//  Created by Yoshi Miyazaki on 2015/04/10.
//  Copyright (c) 2015 Yoshi. All rights reserved.
//

#include "gibbse.h"

extern const int numofElement;


/*--------------------------------------------------------------------------
  class: Molecule_Gibbs
  ---------------------------------------------------------------------------*/
Molecule_Gibbs::Molecule_Gibbs(string _s, double Ti, double Pi) : Molecule(_s){
    // Search for Gibbs Energy
    tensor1d<double> vG(0.0,8);
    T = Ti;  P = Pi;                      /* store T&P in this object  */
    gibbsE = data_search(_s, vG);         /* gibbsE is in [J] NOT [kJ] */
    
    // convert Gibbs energy into given condition
    double Pbar = P/1e5;   /* convert to bar */
    if (phase == 0){gibbsE /= (R*T); gibbsE += log(Pbar);} // gas    phase
    if (phase == 1){gibbsE /= (R*T);                     } // liquid phase
    if (phase == 2){gibbsE /= (R*T);                     } // solid  phase
}
Molecule_Gibbs::Molecule_Gibbs(string _s, double Ti, double Pi, tensor1d<double> _v) : Molecule(_s){    
    // Search for Gibbs Energy
    T = Ti;  P = Pi;                      /* store T&P in this object  */
    gibbsE = data_search(_s, _v);         /* gibbsE is in [J] NOT [kJ] */
    
    // convert Gibbs energy into given condition
    double Pbar = P/1e5;   /* convert to bar */
    if (phase == 0){gibbsE /= (R*T); gibbsE += log(Pbar);} // gas    phase
    if (phase == 1){gibbsE /= (R*T);                     } // liquid phase
    if (phase == 2){gibbsE /= (R*T);                     } // solid  phase
}
double Molecule_Gibbs::data_search(string s_mol, tensor1d<double> _v){
    /* read Gibbs data from various data format */
    string file1  = "/Users/yoshi/Documents/OneDrive/OneDrive/Programs/C++/GibbsFE_minimization/Gibbs_data/Gibbs_Stixrude.txt";  /* first from Stixrude */
    
    string s_name; char c_phase;
    double Gf = 0.0;
    
    /* first check Stixrude database */
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
    
    /* for SiO2(l) use Stixrude's formulation */
    if (s_mol == "SiO2(l)"){
        K0 = _v[0], K0d = _v[1];
        double alpha0 = _v[2], C0 = _v[3], C1 = _v[4], d0 = _v[5];
        Gf = GibbsE_SiO2(K0, K0d, alpha0, C0, C1, d0);
	
        set_phase('l');
        return Gf;
    }    
    
    /* if not found, look for liquid-database file */    
    string file2 = "/Users/yoshi/Documents/OneDrive/OneDrive/Programs/C++/GibbsFE_minimization/Gibbs_data/melting-curve/" + s_mol + ".txt";
    ifstream fmelt(file2);
    if (!fmelt){ cout << "No melting curve data for " << file2 << ". " << endl; exit(8);}
    else{
        Gf = GibbsE_melting(fmelt, _v[0]);
        c_phase = 'l';  set_phase(c_phase);
        fmelt.close();
        return Gf;
    }
    
    return Gf;
}
void Molecule_Gibbs::set_phase(char c_phase){
    if(c_phase=='s'){phase = 2;}
    if(c_phase=='l'){phase = 1;}
    if(c_phase=='g'){phase = 0;}    
}
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
double Molecule_Gibbs::GibbsE_SiO2(double K0, double K0d, double alpha0,
                                   double C0, double C1, double d0){
    /* reference Gibbs Energy: qtz at (Tmo, Patm) */
    double Tm0 = 1710,  Patm = 1.0e5;
    tensor1d<double> v_tmp(0.0,1);
    Molecule_Gibbs   qtz("SiO2(q)",Tm0,Patm,v_tmp);
    double           Gqtz = qtz.getGibbsE()*(R*Tm0);

    /* refernce volume */
    double V0 = 24.5e-6;
    
    /* calculate Gibbs free energy using EoS by Stixrude (2011) */
    double VP = V0*K0/(K0d-1)*(pow(1.0+K0d*P/K0, 1-1/K0d)-1);
    double ST = int_SdT(Tm0, V0, K0, K0d, C0, C1, alpha0, d0);
    
    double Gf = Gqtz + VP - ST;
    
    return Gf;
}
double Molecule_Gibbs::int_SdT(double Tm0, double V0, double K0, double K0d, 
                               double C0,  double C1, double alpha0, double d0){
    /* integrate S from T0 to T */
    /* reference entropy at T=T0, P=Patm */
    double S0 = 47.928;
    double SdT = (S0 - int_alVdP(V0, K0, K0d, alpha0, d0, T))*(T-Tm0);
    if (1){
        SdT += C0*(T*log(T/Tref) - Tm0*log(Tm0/Tref) - (T-Tm0));
        SdT += 0.5*C1*((T-Tref)*(T-Tref) - (Tm0-Tref)*(Tm0-Tref));
        /* C0*log(Td) - C0*log(T0) + C1*(Td-T0)
           reorganzie => 
           C0*(Tlog(T/T0) - T) + 0.5*C1*(T-T0)^2 */
        
    }else{
        SdT += C0*(T*log(T/Tref) - Tm0*log(Tm0/Tref) - (T-Tm0));
        SdT += C1*int_sqrtCpdT(Tm0, T);
    }

    return SdT;
}
double Molecule_Gibbs::int_sqrtCpdT(double Tbeg, double Tend){
    /* integrate SdT where TdS ~ C1*sqrt(T-T0) */
    int N = 200;  double dT = (Tend-Tbeg)/double(N-1);
    
    double SdT = 0, Sa, Sb, Ta, Tb;
    for (int i=0; i<N; i++){
        Ta = Tbeg + dT*(double)(i);
        Tb = Tbeg + dT*(double)(minv(i+1,N-1));
        Sa = 2*sqrt((Ta-Tref)/(Ta*Ta))*Ta*(sqrt(Ta-Tref)-sqrt(Tref)*atan(sqrt(Ta/Tref-1)))/sqrt(Ta-Tref);
        Sb = 2*sqrt((Tb-Tref)/(Tb*Tb))*Tb*(sqrt(Tb-Tref)-sqrt(Tref)*atan(sqrt(Tb/Tref-1)))/sqrt(Tb-Tref);
        SdT += (Sa+Sb)*dT * .5;
        //cout << Sa << " , " << Sb << "  @" << Ta << "   =  " << SdT << endl;
    }
    
    return SdT;
}
double Molecule_Gibbs::int_alVdP(double V0,     double K0, double K0d, 
                                 double alpha0, double d0, double T1){
    /* integrate dS/dP = -dV/dT = -alpha*V from P0 to P */
    int N = 200;    double dP = (P-Patm)/double(N-1);
    
    double dS = 0, Xa, Xb;
    for (int i=0; i<N; i++){
        double Pa = Patm + dP*(double)i;
        double Pb = Patm + dP*(double)(i+1);

        if (Anderson == 1){
            Xa = exp(-1.*d0/K0*Pa)/pow(1+K0d*Pa/K0,1/K0d)*exp(alpha0*(T1-Tref));
            Xb = exp(-1.*d0/K0*Pb)/pow(1+K0d*Pb/K0,1/K0d)*exp(alpha0*(T1-Tref));
            dS += alpha0*V0*(Xa+Xb)*dP/2.0;
        }else{
            Xa = 1./pow(1+K0d*Pa/K0,1/K0d)*exp(alpha0*(T1-Tref))/(1 +K0d*Pa/K0);
            Xb = 1./pow(1+K0d*Pb/K0,1/K0d)*exp(alpha0*(T1-Tref))/(1 +K0d*Pa/K0);
            dS += alpha0*V0*(Xa+Xb)*dP/2.0;
        }
    }
    
    return dS;
}
double Molecule_Gibbs::GibbsE_melting(ifstream& fmelt, double dSP){
    tensor1d<double> v_tmp(0.0,1); /* to create solid Molecule_Gibbs */
    /* read melting temperature from fmelt, and convert into Gibbs free energy data */
    double Gf;
    string sline; /* store each line */
    
    double V0, DS0, A;   /* read the 1st line, molar volume, entropy change at 0Pa,
                            d(dS0)/dP are stored */
    fmelt >> V0 >> DS0 >> A;
    V0 *= 1.0e-6;    /* cm3/mol to m3/mol */
    
    double Tc, Tb = 0, Ta = 0;
    double Pc, Pb = 0, Pa = 0;
    string s_prev = "A", s_solid;
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
    
    /* when T is over the data range,
       extrapolate using last two lines  */
    if (P > Pc){ Tb = Ta; Pb = Pa; } 
    
    /* calculate melting tempeature at T=T from (Tprev, Pprev) and (Tm,Pm) */
    double Tmelt = internal_division(P, Pb, Pc, Tb, Tc);
    // cout << "* - calc melting of " <<  s_prev << endl;
    // cout << "Melting T at: " << setprecision(5) << Tmelt << " K at P = " << P << "[K]. b/w " << Tb << " and " << Tc << " K, P =" << Pb/1e9 << " and " << Pc/1e9 << "GPa. " << endl;
    
    /* Gibbs free energy should be the same b/w solid and liquid at melting temperature.
       Gibbs energy of liquid would be calculated by adding 
       the difference of Gibbs energy from the melting temperature */
    Molecule_Gibbs mole_solid(s_prev, T, P, v_tmp);
    
    /* calc liquid entropy as Sliq = Ssol + dS
       ... use the fact that dS does not change much w/ temperature.
       we assume dS as a simple function of pressure
       dS = dS0 [1 - A*(1 - exp(-P/P0))] 
    
       ... in dataset, dS0 and A must be specified. 
       P0 is set to 100 GPa in our model. */
    double P0 = 50e9;
    double DS = DS0*(1 - dSP*(1 - exp(-1.*P/P0)));
    
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
    // cout << "   G of sol: " << Gsol_T/(R*T) << " - liq: " << Gf/(R*T) << endl<<endl;
    
    return Gf;
}
double Molecule_Gibbs::GibbsE_solid(string s_prev, double Tmelt, double Pmelt){
    tensor1d<double> v_tmp(0.0,1);
    /* calculate GFE of s_prev (solid) at T=Tmelt, P=Pmelt */
    Molecule_Gibbs mole_solid(s_prev, Tmelt, Pmelt, v_tmp);
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


/*--------------------------------------------------------------------------
  class: MoleculeData_G
  ---------------------------------------------------------------------------*/
MoleculeData_G::MoleculeData_G(string filename, double T, double P){
    tensor1d<double> vG(0.0,8);
    vG[0]=0.964;    vG[1]=-0.2861; vG[2]=2.626e+10; vG[3]=4.1252;
    vG[4]=1.549e-5; vG[5]= 3.32;   vG[6]=0.055921;  vG[7]=0.23;
    
    *this = MoleculeData_G(filename, T, P, vG);
    //numofMolecule = database.size();
    //list_species  = filename;
}
MoleculeData_G::MoleculeData_G(string filename, double T, double P, tensor1d<double> v_dG){
    // Open File
    ifstream fMo(filename);
    if (!fMo){ cout << "Error in MoleculeData_G: Failed to open molecule data: " << filename << endl;}
    else{ } // cout << "File Opened (Molecule Data : " << filename << " )" << endl;}
    
    // Read line one by one
    while(!fMo.eof()){
        string s_name;
        fMo >> s_name;
        if (s_name == ""){break;}
        
        /* store dG */
        tensor1d<double> _vp(0.0,6);
        if      (s_name == "MgO(l)"){ _vp[0] = v_dG[0]; }
        else if (s_name == "FeO(l)"){ _vp[0] = v_dG[1]; }
        else if (s_name =="SiO2(l)"){
            for (int i=0; i<6; i++){  _vp[i] = v_dG[i+2]; }}
        Molecule_Gibbs data_molecule(s_name, T, P, _vp);
        database.push_back(data_molecule);
    }
    fMo.close();
    
    numofMolecule = database.size();     // # of molecule
    list_species  = filename;
}
MoleculeData_G& MoleculeData_G::operator=(const MoleculeData_G& copy){
    numofMolecule = copy.numofMolecule;
    list_species  = copy.list_species;
    database      = copy.database;
    molarmass     = copy.molarmass;
    
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
void MoleculeData_G::rm_nonexisting(tensor1d<double>& ninit){
    /* remove species that won't appear in the given set of elemental composition
       e.g.) FeO should be removed if there are no Fe in the system... 
    
       so, check for each considering species (int j, j'th molecule[p] != 0)
       whether it includes any non-existing elements (int p, element_mole[p] == 0) */

    /* element composition ... add up each molecule information */
    tensor1d<double> element_mole(0.0,numofElement);
    for (int p=0; p<numofElement; p++){
        for (int i=0; i<(int)numofMolecule; i++){
            element_mole[p] += ninit[i] * database[i].getComposition(p);
        }
    }
    
    /* decide on the species that should be removed. */
    tensor1d<int> rm_species(0,numofElement);
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
            ninit.erase(j);                     /* adjust n_init      */
        }                   
    }
    // cout << "no. of Molecule = " << numofMolecule << endl;
    // for (int j=0; j<numofMolecule; j++){ cout << database[j].getMoleculeName() << endl; }
}
tensor1d<int> MoleculeData_G::getphase(){
    tensor1d<int> phase(0, numofMolecule);
    for (int j=0; j<(int)phase.size(); j++){ phase[j] = database[j].getphase(); }
    return phase;
}
