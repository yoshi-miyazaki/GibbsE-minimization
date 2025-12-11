/*
  min.cpp
  - calc melt fraction for a given T & P
 */

#include "./CGgibbsmin.h"

int main(){

    /* set T & P through the input stream */
    double T = 2500, P = 30e9;
    cin >> T >> P;

    /* prepare chemical potential for species listed in molcule_simple */
    string           list_txt = "./molecule_simple.txt";
    MoleculeData_G   slist(list_txt, T, P);
    tensor1d<double> n_init(0.0, (int)slist.getnumofMolecule());

    /* set initial composition. 
       ... this model is MgO-FeO-SiO2 ternary system. Input them in molar amount. Not wt% */
    int m_comp = (int)n_init.size();
    n_init[m_comp-3] = 540;   // MgO
    n_init[m_comp-2] = 60;    // FeO
    n_init[m_comp-1] = 500.;  // SiO2
    
    /* results  */
    cout << "Equilibrium compositions at" << endl;
    gibbsminCG        res(slist, n_init);
    tensor1d<double>  nbest = res.getnbest();
    
    cout << fixed << T << " K. " << P/1e9 << " GPa." << endl;
    for (int j=0; j<m_comp; j++){
	cout << slist[j].getMoleculeName() << "\t" << nbest[j] << " \t <- " << n_init[j] << "\t" << slist[j].getGibbsE() << endl;
    }
    cout << endl;

    /* for density calculation */
    double eP = 1e5;
    MoleculeData_G    slistP(list_txt, T, P+eP);
    gibbsminCG        resP(slistP, nbest);
    tensor1d<double>  nbestP = resP.getnbest();
    for (int j=0; j<m_comp; j++){
	cout << slist[j].getMoleculeName() << "\t" << nbestP[j] << " \t <- " << n_init[j] << "\t" << slist[j].getGibbsE() << endl;
    }
    
    double GP = resP.getGbest()*(R*T);
    double G0 =  res.getGbest()*(R*T);
    
    tensor1d<double> mass_spec(0.0,(int)slist.getnumofMolecule());
    for (int i=0; i<m_comp; i++){ mass_spec[i] = slist[i].getWeight(); }
    double V0 = (GP-G0)/eP, umass = n_init*mass_spec, rho = umass/V0;
    cout << "rho = " << rho << endl;
    
    /*
    ifstream fin("./Pgrid200.txt", ios::in);
    ofstream fout("./mf.txt", ios::out);
    int M = 100;
    for (int i=0; i<30; i++){
	double P = (double)i*1e9;

	fin >> P;
	fout << P/1e9 << "\t";
	
	for (int j=0; j<M; j++){
	    double T = 1500 + 4000*j/(M-1);
	    
	    MoleculeData_G  slist(list_txt, T, P);
	    gibbsminCG        res(slist, n_init);
	    tensor1d<double>  nbest = res.getnbest();

	    double mf = res.melt_mfrac();

	    fout << mf << "\t";
	    
	}
	fout << endl;
    }

    fin.close();
    fout.close();
    */
    
    return 0;
}
