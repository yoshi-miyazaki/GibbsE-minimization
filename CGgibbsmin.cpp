//
//  gibbsminCG.cpp
//  Gibbs new version using projection
//
//  Created by Yoshi Miyazaki on 2015/08/04.
//  Copyright (c) 2015 Yoshi. All rights reserved.
//

#include "CGgibbsmin.h"

gibbsminCG::gibbsminCG(MoleculeData_G& moleculelist, tensor1d<double>& n_starting){
    // cout << endl << "--- [Conjugate Gradient] ---" << endl;
    /* store basic information */
    n_in = n_starting;
    list = moleculelist;   list.rm_nonexisting(n_in);
    if (list.getnumofMolecule() == 30){ exit(0); }
    
    /* thermodynamic data */
    gibbse = getGibbslist(list);      /* Gibbs energy (converted to c_j)        */
    phase  = getphaselist(list);      /* phase                                  */
    //cout << "m.size(): " << moleculelist.getnumofMolecule() << " ad: " << &moleculelist << endl;
    
    /* pressure, temperature, and starting vector for conjugate gradient */
    T = list[0].getT();   P = list[0].getP();
    
    // for (int i=0; i<gibbse.size(); i++){cout << " mu[" << i << "] = " << gibbse[i] << endl;}
    
    /* non-ideal melt and solid-solution */
    tensor1d<double> W0(0.0,6);
    if (is_oxide){
        W0[0] = -125e3;     W0[1] = 1.15e-6;   W0[2] = 1;
        W0[3] = -124e3;     W0[4] = 0.34e-6;   W0[5] = 12.8;
        meltSi.add_endmember("MgO(l)","FeO(l)","SiO2(l)");
        meltSi.add_Margules(T,P,W0);
    }
    create_sslist();
    
    /* mass-balance matrix */
    massBalance  massbo(list);
    massm  = massbo.get_massBalance();
    atomic = massbo.get_atomic();
    
    /* conjugate gradient */
    conjugateG();
    //cout << (double)(end-begin)/CLOCKS_PER_SEC << endl;
    
    /* put back the molecule */
    if (moleculelist.getnumofMolecule() != list.getnumofMolecule()){
        tensor1d<double> nbest_resize(0.0, moleculelist.getnumofMolecule());
        tensor1d<int>    phase_resize(0.0, moleculelist.getnumofMolecule());
        for (int j=0; j<list.getnumofMolecule(); j++){
            cout << list[j].getMoleculeName() << "\t" << list.getnumofMolecule() << endl;
            for (int i=0; i<moleculelist.getnumofMolecule(); i++){
                if (moleculelist[i].getMoleculeName() == list[j].getMoleculeName()){
                    nbest_resize[i] = nbest[j];
                    phase_resize[i] = phase[j]; break;
                }
            }
        }
        nbest = nbest_resize;
        phase = phase_resize;
        // cout << "CG: " << nbest.size() << "\t" << phase.size() << "\t list: " << list.getnumofMolecule() << endl;
    }

    /* check element conservation */
    tensor1d<double> qin  = massm.transpose()*n_starting;
    tensor1d<double> qout = massm.transpose()*nbest;
    for (int k=0; k<(int)qin.size(); k++){
        if (abs(qin[k] - qout[k]) > 1){
            nbest = n_starting;
            break;
        }
    }

}
void gibbsminCG::conjugateG(){
    int m = (int)n_in.size();    /* size of system */
    srand((unsigned int)time(NULL));
    
    /*
      cout << "mass balance" << endl;
      int q = massm.nrows(), r = massm.mcols();
      for (int i=0; i<q; i++){
      for (int j=0; j<r; j++){
      cout << setw(5) << setprecision(0) << massm[i][j] << "  ";
      }
      cout << endl;
      }
      cout << "gibbse" << endl;
      int sf = gibbse.size();
      for (int i=0; i<sf; i++){
      cout << setw(20) << list[i].getMoleculeName() << " - " << fixed << setprecision(7) << gibbse[i] << endl;
      }*/

    /* define variables */
    tensor1d<double>  n_old = n_in, n = n_in;
    tensor1d<double> dn_old(0.0,m), dn(0.0,m),  cg(0.0,m);
    double    oldG = gibbsE(n_old), G = oldG,   beta;
    
    /* calculate projection matrix to save computing time. 
       ... when exist = 1 for all i, use projv.
       when exist = 0 for any i, recalculate projection matrix
       and save it in exist_old & exist_saved. */
    projv  = compose_projection(massm);
    projv_saved.resize(0.0,m,m);
    tensor1d<int>  exist(1,m);
    tensor1d<int>  exist_old = exist, exist_saved = exist;
    
    /* convergence parameters */
    int comp = 1;    /* when comp = 0 & count > 700, use zerogas() to calculate grad.
                        After dn.norm() becomes small enough with zerogas(),
                        set comp = 1, and calculate grad in the complete dimension.	    
                        ==> Convergence criteria is dn.norm() < eps (and comp > 1) */
    int sd = 1, up = 0;  /* sd = 1: cg = 0 (steepest descent) in the next step.
                            up = 1: G - oldG > 0.
                            ==> used inside bisection() */
    
    /* termination conditions */
    nbest = n_in;  /* optimized composition is stored in nbest */
    int  count = 0, count_cg = 0, bound = 10, CMAX = 9e2, show = 1e4;
    double eps = 5e-3;
    
    while (count < CMAX){
        if (count > CMAX/2){ eps = 1e-2; }
        
        count++;  count_cg++;
        tensor1d<double> n_save = n_old;
        /* store results from previous step. */
        n_old = n;  dn_old = dn;  oldG = G;  exist_old = exist;
	
        /* delete un-resolvable numbers */
        double n_unres = n_old.max()*(LIM*0.1);
        for (int i=0; i<m; i++){ if(n_old[i] < n_unres && n_old[i] > 0){ n_old[i] = 0.0;}}
        for (int i=0; i<m; i++){ if (std::isnan(n_old[i])){ cout << "error, nan in : " << i << " at step " << count <<endl; }}
	
        /* calculate gradient */
        exist = 1;    dn = grad(n_old,0)*(-1.0);
        /* trick to converge faster. if gas/liq mole i = 0 => exist[i] = 0 */
        if (count <= bound){ comp = 1; }
        if (count > bound && comp == 0){ zerogas(n_old,dn,exist); }
        avoidneg(n_old,dn,exist);  /* project to balance mass & delete unnessary spec */

        /* avoid nan in dn (grad)
           ... `nan` can appear in dn (grad) after operating `mass_balance` 
           .   The program stops the calculation in such case because further calculation is not possible.
        
           .   This happens when one of the element disappears from the system.
           .   Projection matrix used in `mass_balance` will not be calculable, and `nan` appears.
        */
        bool isnan_grad = dn.isnan();
        if (isnan_grad){ converge = false; break; }
        
        /* conjugate gradient
           ... Calculate `beta` and cg.
           .   Note that we use steepest descent when the dimension of vector changes.
           .   beta = 0: conguate => steepest descent */
        if (exist_old != exist){ sd = 1;}
        if (sd == 0){ beta = (dn*(dn-dn_old))/square(dn_old.norm()); } /* calc coeff */
        else        { beta = 0; }
        beta = ((beta<0) ? 0.0 : beta);  /* beta can't theoretically be negative */
	
        /* use steepest decent every m step. 
           ... this is likely to converge the gradient faster! */
        if ((count_cg % m) == 0){beta     = 0;}
        if (beta == 0)          {count_cg = 0;}

       	cg  *= beta;  cg += dn;          /* conjugate gradient   */
        sd   = 0;                        /* reset "sd" indicator */
	
        /* avoid numerical error.
           ... if cg is not 0 and cg < n_old*LIM,
           the difference between n_old & n cannot be resolved. */
        bool cg2small = 0;
        for (int j=0; j<m; j++){
            if (0 < cg[j] && cg[j] < n_old[j]*LIM){ 
                cg2small = 1; // cout << "cg: " << cg[j] << "n: " << n_old[j]*LIM << endl;
            }
        }
	
        /* If comp = 0, use complete dimension from next step. */
        if (cg2small || cg.norm() < eps){         /* if cg is too small (cg2small) */
            if (comp == 0){ comp = 1; continue; }}  /* AND comp = 0, reset the search and
                                                       calculate cg in complete dimension */
	
        /* termination criteria => |dn| < eps for M consecutive times,
           but only if dn is calced in complete dimension (which means comp == 1). */
        if (dn.norm() < eps && comp > 0){
            Gbest = G;
            
            /* REMOVE THIS break; WHEN DEALING WITH MELTS
               ... else, run avoidneg for M times to check whether the current n is the actual optimal value */
            break;
            
            int M = m*10;
            for (comp=1; comp<M; comp++){   /* when i = i, comp = (i+1) */
                exist = 1;  dn = grad(n_old,1)*(-1.0);  /* with mixck flag on! (1) */
                avoidneg(n_old,dn,exist);
		
                /* check whether any direction is heading towards smaller G */
                double  dc = 1;   n = shrink(n_old,dn,dc);
                bool small = ck_smaller(n_old,n,dn,exist);
                if (dn.norm() > eps && small){ cg = dn; break; }
            }

            /* If |dn| < eps for M times, the current n is likely to be the optimal. */
            if (comp == M){
                cg = dn;
                // cout << "completed" << endl;
                break;
            }
        }
        comp = 0; /* reset comp here. */
	
        /* 
           check the sign of cg. 
           ... even if n[i] = 0 and dn[i] >=0, cg[i] < 0 may occur.
           .   in such cases, reset beta = 0 => cg = dn
        */
        for (int i=0; i<m; i++){
            if (n_old[i]==0.0 && cg[i]< 0){ beta = 0.0;  cg = dn; break; }}
        
        if (count > show){ /* output to check result */
            cout << endl << count << " ... comp: " << comp;
            cout << ", up: " << up << ", sd?: " << sd << endl;
            cout << "G - oldG : " << fixed << G-oldG << ", beta: " << beta <<  endl;
            for (int i=0; i<m; i++){
                cout << setw(15) << list[i].getMoleculeName() << " -old: " << setw(10) << setprecision(7) << scientific << n_old[i] << " - dn = " << setw(10) << fixed << dn[i] << " - cg: " << cg[i] << " - ex: " << exist[i] << endl;
            }
            cout << "done. " << endl << endl;
	    
            cout << "chemi pote: " << endl;
            tensor1d<double> mup = grad(n_old,0);
            for (int i=0; i<(int)mup.size(); i++){
                cout << i << " - " << mup[i] << endl;
            }
        }
	
        /* next point / modify(*dn) : factor to adjust neg mole   */
        double mod;
        n = shrink(n_old,cg,mod);              /* shrink arrow to avoid neg mole       */
        tensor1d<double> nright = n;
        bisection(n_old,n,cg,mod,exist,sd,up); /* bisection or golden section search   */
        nbest  = n;  G = gibbsE(nbest);        /* when dG = 0 is not found, sd will be
                                                  updated to 1, and the beta next step 
                                                  is set to 0. (sd = steepest descent) */
        
        /*if (count > show){
            tensor1d<double> gc = grad(nbest,0);	
            cout << endl;
            cout << count << " ..c = " << comp << ", up: " << up << ", sd?: " << sd << endl;
            cout << "G - oldG : " << fixed << G-oldG << ", beta: " << beta <<  endl;
            for (int i=0; i<m; i++){
                cout << setw(15) << list[i].getMoleculeName() << " -old: " << setw(10) << setprecision(7) << scientific << n_old[i] << " -new: " << n[i]  << "  - shrinked: " << nright[i] << " - dn = " << setw(10) << fixed << dn[i] << " - cg: " << cg[i] << " = " <<  fixed << exist[i] << " - mu = " << gc[i] << endl;
            }
            cout << "element mass balance" << endl;
            tensor1d<double> qold = massm.transpose()*n_old;
            tensor1d<double>    q = massm.transpose()*n;
            for (int j=0; j<(int)qold.size(); j++){ cout << qold[j] << "  " << q[j] << endl; } 
	    
            cout << "done. " << cg.norm() << endl << endl;
            }*/
    }
    
    /* output the final result. 
       cout << endl << count << "---------------------- |dG/dn| = " << dn*dn << "\t @" << T << endl;
       cout << "- n:          cg:              dn:            exist:" << endl;
       for (int j=0; j<m; j++){
       cout << setw(17) << nbest[j] << "  " << setw(16) << cg[j] << "  " << setw(16) << dn[j] << " " << exist[j] << endl;}
       cout << "absdn : " << dn.norm() << endl << "count: " << count << endl; */
    
    /*
      end of conjugateG()
      ... store the result of convergence to `converge`
    */
    converge = (count == CMAX) ? false : true;
}
tensor1d<double> gibbsminCG::getGibbslist(MoleculeData_G& moleculelist){
    int m = moleculelist.getnumofMolecule();
    
    tensor1d<double> glist(m);
    for (int i=0; i<m; i++){ glist[i] = moleculelist[i].getGibbsE(); }
    return glist;
}
tensor1d<int>    gibbsminCG::getphaselist(MoleculeData_G& moleculelist){
    int m = moleculelist.getnumofMolecule();
    
    tensor1d<int> plist(m);
    for (int i=0; i<m; i++){ plist[i] = moleculelist[i].getphase(); }
    return plist;
}
double gibbsminCG::gibbsE(tensor1d<double>& n){
    int m = (int)n_in.size();            /* size of system               */
    int n_mixing = (int)ss_list.size();  /* no. of solid solution models */
    /* delete all mole & configuration entropy in each "solid_solution" obejct */
    for (int k=0; k<n_mixing; k++){ ss_list[k].reset(); }
    
    /* calculate total mole of gas&liquid + store amount in "solid_solution" object */
    double Ngas = 0.0, Nliq = 0.0;
    for (int j=0; j<m; j++){
        if      (phase[j] == 0){ Ngas += n[j]; }
        else if (phase[j] == 1){ Nliq += n[j];
            string _s = list[j].getMoleculeName();
            meltSi.add_mole(_s, n[j]);
        } else {
            string _s = list[j].getMoleculeName();
            for (int k=0; k<n_mixing; k++){
                if (ss_matrix[j][k]){ ss_list[k].add_mole(_s,n[j]); }
            }
        }
    }
    
    double G = 0.0;         /* calculate Gibbs free energy */
    for (int j=0; j<m; j++){
        double dG = n[j]*gibbse[j]; /* if usual solid, simply add mole*mu0 */
        if (n[j] > 0){
            if      (phase[j] == 0){ dG = n[j]*(gibbse[j] +log(n[j]) - log(Ngas));}
            else if (phase[j] == 1){
                string _s = list[j].getMoleculeName();
                double dGmix = meltSi.activity(_s);
                dG = n[j]*(gibbse[j] +log(n[j]) - log(Nliq) + dGmix);
            }
            else {
                string _s = list[j].getMoleculeName();
                for (int k=0; k<n_mixing; k++){
                    if (ss_matrix[j][k]){ dG+= n[j]* ss_list[k].config(_s, 0.0);}
                }/* end of solid solution model search */
            }/* end of gas/liquid/ss search */
        }
        G += dG;  /* add to total Gibbs free energy */
    }
        
    return G;
}
tensor1d<double> gibbsminCG::grad(tensor1d<double>& n, int mixck){
    int m = (int)n_in.size();            /* size of system               */
    int n_mixing = (int)ss_list.size();  /* no. of solid solution models */
    /* delete all mole & configuration entropy in each "solid_solution" obejct */
    for (int k=0; k<n_mixing; k++){ ss_list[k].reset(); }    
    
    /* calculate mixing of gas, liquid, and solid-solution (ss) */
    double Ngas = 0.0, Nliq = 0.0;
    tensor1d<double> Nss(0.0, n_mixing); /* store total mole of each ss system */
    for (int j=0; j<m; j++){
        string _s = list[j].getMoleculeName();
        if      (phase[j] == 0){ Ngas += n[j]; }
        else if (phase[j] == 1){ Nliq += n[j];
            meltSi.add_mole(_s, n[j]);
        }else{
            for (int k=0; k<n_mixing; k++){
                if (ss_matrix[j][k]){ ss_list[k].add_mole(_s,n[j]); }
            }
        }
    }
    
    /* calculate dG gradient in a space*/
    tensor1d<double> dG(m);
    double           non0 = abs(n.max()) * LIM;
    for (int j=0; j<m; j++){
        string _s = list[j].getMoleculeName();
        dG[j] = gibbse[j];
        /* for gas, liquid, and solid-solution, assume ideal mixing
           ... chemical potential should be modified as mu = mu0 + RT*log(concentration)
           BUT for very small concentration, smaller than the numerical digit 
           than C++ could handle, assume non0 mole instead of actual amount.
           This would prevent large gradient value and make converge faster */
        if (phase[j] == 0){
            if (n[j] > 0)    { dG[j] = gibbse[j] + log(n[j]) - log(Ngas);}
            else if(Ngas > 0){ dG[j] = gibbse[j] + log(non0) - log(Ngas);}
        }else if(phase[j] == 1){
            double dGmix = meltSi.activity(_s);
            if (n[j] > 0)    { dG[j] = gibbse[j] + log(n[j]) - log(Nliq) + dGmix;}
            else if(Nliq > 0){ dG[j] = gibbse[j] + log(non0) - log(Nliq) + dGmix;}
        }else {
            dG[j] = gibbse[j];
            for (int k=0; k<n_mixing; k++){
                if (ss_matrix[j][k]){
                    double Ntot  = ss_list[k].tot_mole();
                    //cout << " -- " << k << " , " << ss_list[k].config(_s) << endl;
                    if     (n[j] > 0){ dG[j] += ss_list[k].config(_s, 0.0); }
                    else if(Ntot > 0){ dG[j] += ss_list[k].config(_s,non0); }
                    // cout<<"n0: "<< _s << " = " << scientific <<  ss_list[k].config(_s,non0) << " , non0 = " << non0 << " , nmax = " << n.max() << endl; }
                    break;
                }
            }
        }
    }
    
    if (mixck == 1){
        if (0){
            // Nliq == 0){
            /* for gas, liquid, and solid-solution,
               when their total amount is 0, no lowering through activity happens. 
               A different phase may be showing a lower chemical potential, but
               these phases may be more stable when mixing is included.
               
               when (ssck == 1), the program calculate Mg# and use the estimated 
               liquid concentration for liquid acitivity. */
            Matrix<double>    trans = massm.transpose();
            tensor1d<double>  q_elm = trans*n;
            
            /* obtain Mg, Fe, Si element mole */
            int eMg = -1, eSi = -1, eFe = -1;
            for (int l=0; l<(int)q_elm.size(); l++){
                if (atomic[l] == 12){ eMg = l; }
                if (atomic[l] == 14){ eSi = l; }
                if (atomic[l] == 26){ eFe = l; }
            }
    
            /* create arbitrary solid solution using rand()
               ... Si(l) only emerges through decomposition, therefore 
               calculated based on mass balance */
            /* use the fact that Mg is more incompatible than Fe. 
               rMg (Mg" of liquid) should be LOWER than the whole system Mg# */
            double eratio_Mg = q_elm[eMg]/q_elm.sum(), rMg = 1, rFe = 1, rSi = 0;
            int iMgOp = list.intspec("MgO(p)"), iMgOl = list.intspec("MgO(l)");
            int iFeOp = list.intspec("FeO(p)"), iFeOl = list.intspec("FeO(l)");
	    
            /* randomize Mg concentration in liquid... if both MgO and FeO are present */
            // srand((unsigned int)time(NULL));
            if (iMgOl > 1 && iFeOl > 1){
                double Dmu  = gibbse[iMgOp]-gibbse[iMgOl],  emin_Mg = eratio_Mg*exp(Dmu);
                if (Dmu < 0){
                    while (rMg < emin_Mg || eratio_Mg < rMg){
                        rMg = (double)rand()/RAND_MAX; rFe = 1 - rMg; }
                } else{ rMg = eratio_Mg; rFe = 1 - rMg; }
                // cout << "emin_Mg: " << emin_Mg << " - eratio_Mg: " << eratio_Mg << endl;
                /* ... if liquid MgO has lower chemical potential than the MgO-periclase,
                   solid is clearly unstable in this situation.
                   simply assign the system Mg# as Mg concentration in liquid */
            }
            if (eSi > 0){ rSi = 1.;  /* when Si exists... */
                while (rSi > (q_elm[eSi]/q_elm.sum())){ rSi = (double)rand()/RAND_MAX;}}
	    
            for (int j=0; j<m; j++){
                dG[j] = 0;  string _s = list[j].getMoleculeName();
                // cout << "mixing effect for " << _s << endl;
		
                double A = -1.*n.absnon0min();
                if (_s=="MgO(l)" ){ //int _i = list.intspec("MgO(p)");
                    dG[j] = A*rMg *(1-rSi);}
                if (_s=="FeO(l)" ){ //int _i = list.intspec("FeO(p)");
                    dG[j] = A*rFe *(1-rSi);}
                if (_s=="SiO2(l)"){ dG[j] = A*rSi; } // gibbse[j] + log(q_elm[eSi]/q_elm.sum()); }
            }
            // cout << "rMg: " << rMg << " , rFe: " << rFe << " , rSi: " << rSi << endl;

            if (std::isnan(dG[2]) || std::isnan(dG[3])){
                cout << "MgO(l): " << dG[2] << endl;
                cout << "FeO(l): " << dG[3] << endl;
                cout << "rMg = " << rMg << " , rFe = " << rFe << " , rSi = " << rSi << endl;
            }
        }
    }
    return dG;
}
void gibbsminCG::mass_balance(tensor1d<double>& dn, tensor1d<int>& exist){
    /* size of mass-balance matirx */
    int m = massm.nrows();
    int e = massm.mcols();
    tensor1d<int> ones(1,m);

    if (ones == exist){ dn = projv*dn; }
    else if (exist_saved != exist){   /* use same projv if diminish is the same as prev.*/
        Matrix<double> modifiedmassm = massm;
        for (int i=0; i<m; i++){
            if (exist[i] == 0){
                for (int j=0; j<e; j++){ modifiedmassm[i][j] = 0.0;}}
        }
        Matrix<double> modifiedprojv = compose_projection(modifiedmassm);
        dn          = modifiedprojv*dn;
        projv_saved = modifiedprojv;     exist_saved = exist;
    } else {     dn = projv_saved*dn; }
    
    /* if (std::isnan(dn[0])){
       int r = projv.nrows(), q = projv.mcols();
       cout <<"projection" << endl;
       for (int i=0; i<r; i++){
       for (int j=0; j<q; j++){ cout << setw(10) << setprecision(5) << projv[i][j] << "  "; }
       cout << endl;
       }} */
    
    /* assure dn[j] = 0 for diminish[j] = 0 */
    for (int j=0; j<m; j++){
        if (exist[j] == 0){ dn[j] = 0; }
    }
}
Matrix<double> gibbsminCG::compose_projection(Matrix<double>& mass){
    /* size of mass-balance matrix */
    int m = massm.nrows();
    
    Matrix<double> trans = mass.transpose();
    Matrix<double> BTB   = trans*mass;
    Matrix<double> luBTB = BTB.lu_decomp();
    Matrix<double> inBTB = luBTB.lu_inverse();
    Matrix<double> BBTBT = (mass*inBTB)*trans;
    
    Matrix<double> proj(0.0,m,m);
    for (int i=0; i<m; i++){ proj[i][i] = 1.0; }
    proj -= BBTBT;
    proj.numeric0(LIM);
    
    return proj;
}
void gibbsminCG::zerogas(tensor1d<double>& n_old, tensor1d<double>& dn, tensor1d<int>& exist){
    int m = (int)n_in.size();    /* size of system */
    int n_mixing = (int)ss_list.size();
    /* delete all mole & configuration entropy in each "solid_solution" obejct */
    for (int k=0; k<n_mixing; k++){ ss_list[k].reset(); }
    
    /* calculate mixing of gas, liquid, and solid-solution (ss) */
    double Ngas = 0.0, Nliq = 0.0;
    tensor1d<double> Nss(0.0, n_mixing); /* store total mole of each ss system */
    for (int j=0; j<m; j++){
        if      (phase[j] == 0){ Ngas += n_old[j]; }
        else if (phase[j] == 1){ Nliq += n_old[j]; }
        else {
            string _s = list[j].getMoleculeName();
            for (int k=0; k<n_mixing; k++){
                if (ss_matrix[j][k]){ ss_list[k].add_mole(_s,n_old[j]); }}
        }
    }
    
    double non0 = abs(n_old.max()) * LIM;
    for (int j=0; j<m; j++){
        if(n_old[j] == 0){
            if      (phase[j]==0){ dn[j]=0.0; exist[j]=0; }
            else if (phase[j]==1){ dn[j]=0.0; exist[j]=0; }
            else { string _s = list[j].getMoleculeName();
                for (int k=0; k<n_mixing; k++){
                    if (ss_matrix[j][k]){
                        double Ntot = ss_list[k].tot_mole();
                        if (Ntot > 0){dn[j]=0.0; exist[j]=0; break; }
                    }
                }/* end of solid-solution model search */
            }
        }
    }
    
    /* At least one of the component in exist has to be 1 for each element, 
       else the mass_balance matrix is going to be singular */
    for (int k=0; k<massm.mcols(); k++){
        int in = 0, iel = -1;
        for (int j=0; j<m; j++){
            if(massm[j][k] != 0){  /* if species j includes element k, */
                if ( exist[j] != 0 ){ in++; break; } /* and exist[j] = 1, no problem        */
                else { iel = j; }                    /* but exist[j] = 0 for all j, no good */
            }
        }
        if (in == 0){ exist[iel] = 1; }
    }
}
void gibbsminCG::avoidneg(tensor1d<double>& n_old, tensor1d<double>& dn, tensor1d<int>& exist){
    int m = (int)n_in.size();    /* size of system */
    
    tensor1d<double> dt = dn;
    int count = 0, ifneg = -2;
    while (ifneg < 0){
        dt = dn;  count++;
        for (int j=0; j<m; j++){ if(exist[j] == 0){dt[j] = 0.0;} }
        mass_balance(dt,exist);         /* project to meet mass balance eq.  */	
        
        /*for (int i=0; i<m; i++){ if (isnan(dt[i])){
          cout << "--- after: " << count << endl;
          for (int j=0; j<m; j++){
          cout << setw(15) << list[j].getMoleculeName() << " - " << setw(10) << setprecision(10) << fixed << n_old[j] << "  -  " << setw(20) << dt[j] << " - " << exist[j] << endl;
          }
          break;}}*/
        
        /* check whether it satisfies dt[j] >= 0 when n_old[j] = 0.0 */
        ifneg = 1;
        int s = rand()%m;
        for (int j=0; j<m; j++){
            if (n_old[s] == 0 && dt[s] < 0){ ifneg = -1; exist[s] = 0; break;}
            s++; s=s%m;
        }
    }
    dn = dt;
}
tensor1d<double> gibbsminCG::shrink(tensor1d<double>& n_old, tensor1d<double>& dn, double& mod){
    /* adjust dn so that n = n_old+dn have non-negative amount */
    int m = (int)n_in.size();    /* size of system */
    tensor1d<double> n  = n_old +dn;
    tensor1d<double> dt = dn;
    
    /* shrink arrow */
    int dore = -1; double temp;
    mod = 1.0;
    for (int j=0; j<m; j++){
        if (n_old[j] > 0 && n[j] < 0){
            temp = abs(n_old[j])/(abs(n_old[j])-n[j]);
            if (temp < mod){mod = temp; dore = j;}
        }
    }
    dt *= mod;    n = n_old +dt;
    if (dore != -1){ n[dore] = 0.0;}    /* assure n[dore] = 0 */
    
    /* reassure that no negative mole exists */
    double negLIM = abs(n_old.max())*(-1.0*LIM);
    for (int j=0; j<m; j++){
        if (negLIM < n[j] && n[j] < 0){ n[j] = 0.0; }
        // if (n[j] < 0){ cout << "negative mole >_< : " << j << " is " << scientific << n[j] << endl;}
    }
    
    return n;
}
double gibbsminCG::dG_direction(tensor1d<double>& n, tensor1d<double> cg){
    /* calculate change in GFE along cg direction. dG/da = (nabla G)*cg */
    tensor1d<double> n_pert = perturb_zero(n,cg);
    // clock_t end = clock();
    tensor1d<double> nabG   = grad(n_pert,0);
    // clock_t end2 = clock();
    // cout << "A2: " << (double)(end2-end)/CLOCKS_PER_SEC << endl;
    double dGda = std::inner_product(nabG.begin(), nabG.end(), cg.begin(), 0.0);
    if (std::isnan(dGda)){
        cout << "nan... @" << T << " K" << endl;
        for (int j=0; j<(int)n.size(); j++){
            cout << j <<" - "<< n_pert[j] <<" \t cg: "<< cg[j] <<" , grad: "<< nabG[j] << endl;
        }
        
        Matrix<double>    trans = massm.transpose();
        tensor1d<double>  q_elm = trans*n_in;
        for (int k=0; k<(int)q_elm.size(); k++){ cout << q_elm[k] << "\t"; }
        cout << endl;
        tensor1d<double>  q_now = trans*n_pert;
        for (int k=0; k<(int)q_now.size(); k++){ cout << q_now[k] << "\t"; }
        cout << endl;

        exit(2);
    }
    
    return dGda;
}
bool gibbsminCG::ck_smaller(tensor1d<double>& n_old, tensor1d<double>& n,
                            tensor1d<double>& dn, tensor1d<int>& exist){
    /* judge whether direction dn is likely to lower GFE 
       ... check using the change in sign of r = (nabla G)*(dn) */
    bool   Gsmaller = 0;    
    double oldG = gibbsE(n_old), G = gibbsE(n);
    if (oldG > G){ Gsmaller = 1; }
    else { /* if bisection search is possible, 1 */
        double oldr = dG_direction(n_old,dn);
        double    r = dG_direction(n,    dn);
        if (oldr*r < 0){ Gsmaller = 1;}
	
        //double refr = dG_direction(n_old,exist,dn);
    }
    
    return Gsmaller;
}
tensor1d<double> gibbsminCG::perturb_zero(tensor1d<double>& n, tensor1d<double>& dn){
    tensor1d<double> n_pert = n;
    
    double mindn = 0.0; /* perturb zero-mole species to non0 if dn(i) != 0 */
    for (int i=0; i<(int)n.size(); i++){ if (n[i] == 0 && dn[i] > 0){ mindn = max(mindn,dn[i]); }}
    if (mindn == 0.0){ return n_pert; }

    /* reject any negative mole */
    double  non0 = DBL_MIN*1e5;
    n_pert += dn*(non0/mindn);
    // cout << "ratio: " << scientific << non0/mindn << endl;
    
    double mod = 1.0;
    for (int i=0; i<(int)n.size(); i++){
        if (n_pert[i] < 0){ mod = minv(mod, n[i]/(n[i]-n_pert[i]));}
        // cout << "pert: " << (non0/mindn)*mod << " , n-n_pert : " << n[i]-n_pert[i] << endl;
    }
    n_pert = n + dn*(non0/mindn)*mod*(1-DBL_MIN);
    
    return n_pert;
}
void gibbsminCG::bisection(tensor1d<double> n_old, tensor1d<double>& n, tensor1d<double>& cg, double mod, tensor1d<int>& exist, int& sd, int& up){ /* cg: direction (before shrink) */
    int m = (int)n_in.size();     /* size of system  */
    int e = massm.mcols();   /* no. of elements */
    
    /* step-1 : calc r = grad(n)*cg. 
       If r and oldr is different signs, local minimum must exist between (n_old,n)
       
       BUT, when n_gas/n_liq is zero, r may not be calculated correctly.
       for eg, when both nold_(MgO) and nold_(FeO) are 0.0, but dn_(MgO),(FeO) > 0:
       g_old = -grad(n_old) do NOT include the effect of mixing (log(ni/N)),
       but g = -grad(n) DO include,
       and signs may not necessary reflect the existence of local minimum.
      
       In order to avoid that, we perturb n_old a little and define n_pert,
       so that n_pert includes the effect of mixing, yet its value is so close to n_old */
    tensor1d<double> n_out = n;
    tensor1d<double> g = grad(n,0);   mass_balance(g,exist);
    double  oldr = dG_direction(n_old, cg);
    double     r = dG_direction(n,     cg*(-1.))*(-1.);
    
    if(std::isnan(n_old[0])){ cout << "NAN before bisection. n_old(0)= "  << n_old[0] << endl;}
    if(std::isnan(n[0])){ cout << "NAN before bisection. n(0)= "  << n[0] << endl; }
    
    /* step-2 : Decide on the section for line search.
       double (dn) if necessary.
       (dn) may be so small that (n,n+dn) doesn't include local minimum */
    int toosmall = 0;  /* (dn) may violate mass-balance if it's too small compared to (n) */
    if (oldr*r > 0 && mod == 1.0){
        /* qb, qa: molar amount of each element. */
        Matrix<double>   trans = massm.transpose();
        tensor1d<double> q_old = trans*n_old;
        tensor1d<double> q     = trans*n;
	
        /* expand cgs double to avoid searching the same direction
           over and over again. */
        tensor1d<double> cgs = cg;
        while (oldr*r > 0){
            n += cgs; cgs = n-n_old;   /* take new point (n) using cgs. */
	    
            /*confirm that each element is conserving its amount */
            for (int k=0; k<e; k++){
                if(abs(q_old[k]-q[k]) > (q_old[k]*LIM*100)){ n = n_old +cg/2; toosmall=1; }
            }
            if (toosmall == 1){ break; // cout << "too small: elements:" << endl; break; 
                // for (int k=0; k<e; k++){
                // cout << "    --- compare: " << abs(q_old[k]-q[k]) << " - ";
                // cout <<  q_old[k]*LIM << endl; }
            }
	    
            /* break loop if any of the species are negative */
            int ifneg = 1;
            for (int j=0; j<m; j++){ if(n[j] < 0){ ifneg = -1; }}
            if (ifneg < 0){
                n = shrink(n_old, cgs, mod);
                r = dG_direction(n, cg*(-1.))*(-1.);   break;
            }
            
            /* for the next loop */
            r = dG_direction(n, cg*(-1.))*(-1.);  /* renew the value of r */
            
        }
    }
    // cout << " -- pert, r signs: left: " << oldr << " , right: " << r << endl;
    /* (just in case:) check if (n) diverged by doubling (cg) */
    for (int j=0; j<m; j++){
        if (std::isnan(n[j])){ cout << list[j].getMoleculeName() << " has diverged during cg amplification... " << n[j] << endl; exit(10);}}

    /* step-3. actual LINE SEARCH : golden seciton search or biseciton
       ... res = absnon0 min value of vector |n-n_old| */
    tensor1d<double> nleft = n_old, nright = n, ds = n-n_old;
    double res  = 2.0,  oldres = 0.0;
    double oldG = gibbsE(n_old), G = gibbsE(n);
    
    if (oldr*r > 0){ /* the same sign => golden section search */
        //cout << "gss, old= " << oldr << " , r= " << r << " , G-old= " << G-oldG << endl;
        tensor1d<double> nA, nB, g_old, g;
        double GA, GB, tau = 0.61803398875;
        int refc = 0;
        while (oldres != res && refc < 150){
            refc++;
            ds = n - n_old;
            nA = n_old + ds*(1-tau);   GA = gibbsE(nA);
            nB = n_old + ds* tau;      GB = gibbsE(nB);
	    
            if (min(GA,GB) < min(oldG,G)){      /* switch */
                if (GA > GB){n_old= nA; oldG= GA; } /* next search b/w (nA   ,n ) */
                else        {n    = nB;    G= GB; } /*             b/w (n_old,nB) */
            } else { n_out = n; break; }
	    
            /* calculate gradient */
            oldr = dG_direction(n_old,cg);
            r    = dG_direction(n,    cg*(-1.))*(-1.);
            if(oldr*r < 0){ nright = n; cout << "gss: r- " << r << endl; break; }
	    
            /* termination condition*/
            ds = n - n_old;  n_out = (n_old+n)/2.0;
            oldres = res;    res = ds.absnon0min();
        }
        if (oldr*r > 0){ 	/* if minimization fails : next step => sd */
            if(G <= oldG){n_out = n;    sd= 1; up= 0;}
            else{
                if(up==0){n_out= n_old; sd= 1; up= 1;} /* if prev step: cg-> redo w/ sd */
                else     {n_out= n;     sd= 1; up= 0;} /*            sd-> accept G>oldG */
            }
        }
    }
    if (oldr*r < 0){ /* different signs => bisection search */
        // cout << "bisec, old= " << oldr << " , r= " << r << " , G-old= " << G-oldG << endl;
        tensor1d<double> nA, gA;
        double rl = oldr, rr = r, rA;
	
        int refc = 0;  ds = n-n_old;
        while (oldres != res && refc < 100){
            /* refc > 100 should be unnecessary... */
            refc++;
            nA = n_old + ds*0.5;
            rA = dG_direction(nA,cg);
            
            if (rA == 0)  {n = nA;  n_out = n;  break;}
            if (rA*rl < 0){n = nA;  rr = rA;}
            else      {n_old = nA;  rl = rA;}
            
            /* termination condition*/
            ds = n-n_old;
            oldres = res;  res = ds.absnon0min();  /* termination condition */
        }
        n_out = (n_old+n)*0.5;
        up = 0;
        
        /* if rA is not even close to 0 (= when bisection failes) */
        if (abs(rA)>1){
            // cout << " bisection fail oldr = " << rl << " , rA = " << rA << " , rr = " << rr << endl;
            // oldr & rA same sign => real solution is r side
            if (oldr*rA > 0){ n_out = n; }
            else            { n_out = n_old; }
        }
        if (0) { //nleft == n_old || nright == n && 0){
            cout << "match, refc:" << refc << ", oldres: " << scientific << oldres << " res = " << res << " , rA = " << rA << endl;
            tensor1d<double> mcg = cg*(-1.0);
            tensor1d<double> n_pl = perturb_zero(nleft,cg), n_pr = perturb_zero(nright,mcg);
            double r1 = dG_direction(n_pl,cg), r2 = dG_direction(n_pr,cg);
            cout << "rl = " << r1 << " , rr = " << r2 << endl;
            cout << "nFeO, left = " << nleft[3] << " , pl = " << n_pl[3] << " , A = " << nA[3] << " , n_pr = " << n_pr[3] << " , rght = " << nright[3] << endl;

            sd = 1;
            for (int i=0; i<(int)ds.size(); i++){
                cout << " - " << list[i].getMoleculeName() << " = " << ds[i] << " , lef = " << nleft[i] << " n_old: " << n_old[i] << " , lp = " << n_pl[i] << " , A = " << nA[i] << " , rp = " << n_pr[i] << " , right = " << nright[i] << endl;
            }
        }
    }
    n = n_out;
    if(std::isnan(n[0])){ cout << "NAN happened in bisection." << endl; exit(17); }
}
/* output result */
void gibbsminCG::result(double T, double P, string& filename){ /* T in K, P in Pa */
    int m = (int)n_in.size();     /* size of system */
    
    Gbest = gibbsE(nbest);
    ofstream fout(filename, ios::app);
    fout << T << " " << P/1e9 << " " <<setprecision(10) << G;
    for (int j=0; j<m; j++){
        fout << " " << setprecision(10) << nbest[j];
    }
    fout << endl;
    fout.close();
    
    tensor1d<double> mu = grad(nbest,0);
    cout << "------------( Result )------------" << endl << endl;
    cout << "At T = " << T << " K. or " << T - 273.15 << " C ... " << converge << endl;
    cout << "At P = " << P/1e9 << " GPa." << endl;
    cout << "Minimum Gibbs free energy is " << endl;
    cout << "  G = " << fixed << setprecision(12) << Gbest << endl;
    cout << "Composition: " << endl;
    for (int j=0; j<m; j++){
        cout << setw(3) << j << setw(20) << list[j].getMoleculeName() << " phase: " << phase[j] << ",mol = " << fixed << setprecision(10) << nbest[j] << " , mu: = " << mu[j] << endl;
    }
    cout << endl;
}
tensor1d<double> gibbsminCG::getngas_or_solid(int phase_no){
    int m = (int)nbest.size();
    tensor1d<double> gas_or_solid(0.0,m);
    for (int j=0; j<m; j++){
        if (phase[j]==phase_no){ gas_or_solid[j] = nbest[j];}
    }
    if (m<24){ cout << "P2:\t" << gas_or_solid.size() << "\t" << nbest.size() << "\t" << n_in.size() << endl; }
    return gas_or_solid;
}
double gibbsminCG::meltfrac(){
    int m = (int)n_in.size();
    tensor1d<double> mass(0.0,m);
    for (int j=0; j<m; j++){ double _d = list[j].getWeight();  mass[j] = _d;}
    // cout << "name: " << list[j].getMoleculeName() << " - weight: " << _d << endl; 
    
    tensor1d<double> sol(0.0,m), liq(0.0,m), gas(0.0,m);
    for (int j=0; j<m; j++){
        if (phase[j]==2){ sol[j] = nbest[j];}
        if (phase[j]==1){ liq[j] = nbest[j];}
        if (phase[j]==0){ gas[j] = nbest[j];}
    }
    double m_mel = mass*liq, m_sol = mass*sol;        
    return m_mel/(m_mel+m_sol);
}
double gibbsminCG::meltfrac(double T, double P, string& filename){
    int m = (int)n_in.size();
    tensor1d<double> mass(0.0,m);
    for (int j=0; j<m; j++){ double _d = list[j].getWeight();  mass[j] = _d;}
    // cout << "name: " << list[j].getMoleculeName() << " - weight: " << _d << endl; 
    
    tensor1d<double> sol(0.0,m), liq(0.0,m), gas(0.0,m);
    for (int j=0; j<m; j++){
        if (phase[j]==2){ sol[j] = nbest[j];}
        if (phase[j]==1){ liq[j] = nbest[j];}
        if (phase[j]==0){ gas[j] = nbest[j];}
    }
    double m_mel = mass*liq, m_sol = mass*sol;
    
    ofstream fout(filename, ios::app);
    fout << P << " " << T << " " <<setprecision(10) << m_mel/(m_mel+m_sol) ;
    fout << endl;
    fout.close();
    
    return m_mel/(m_mel+m_sol);
}


/*--------------------------------------------------------------------------
// Def of class: solid solution
---------------------------------------------------------------------------*/
void gibbsminCG::create_sslist(){
    /* define the solid-solution system. 
       ... when defining the object, specify <string> of end members, and constitution no.
       e.g.) for olivine, Mg2SiO4(f) will be 2, because of 2 Mg-sites. */
    int n_solution = 14;
    tensor1d<solution>  ssref(n_solution);
    
    solution     olivine("Mg2SiO4(o)", "Fe2SiO4(o)", 2);    ssref[0] = olivine;
    solution  wadsleyite("Mg2SiO4(w)", "Fe2SiO4(w)", 2);    ssref[1] = wadsleyite;
    solution ringwoodite("Mg2SiO4(r)", "Fe2SiO4(r)", 2);    ssref[2] = ringwoodite;
    solution  perovskite("MgSiO3(p)" , "FeSiO3(p)" , 1);    ssref[3] = perovskite;
    solution    mgwusite("MgO(p)"    , "FeO(p)"    , 1);    ssref[4] = mgwusite;
    solution    melilite("Ca2Al2SiO7","Ca2MgSi2O7" , 1);    ssref[5] = melilite;
    solution      spinel("Mg4Al8O16(s)","Fe4Al8O16(s)",4);  ssref[6] = spinel;
    solution     ferrite("MgAl2O4(c)"  ,"FeAl2O4(c)" , 1);  ssref[7] = ferrite;
    solution opxM2("MgSiO3(o)",0.5,1,"FeSiO3(o)",0.5,2,"MgAl2SiO6(o)",1,3,1);
    solution cpxM1("CaMgSi2O6(c)"  ,"CaFeSi2O6(c)"    ,"CaAl2SiO6(c)",    1);
    solution garnetX("Mg3Ai2Si3O12(g)",1,1,"Fe3Al2Si3O12(g)",1,2,"Ca3Al2Si3O12(g)",1,3,"Mg4Si4O12(g)",1,1,3);
    solution garnetY("Mg3Ai2Si3O12(g)",1,1,"Fe3Al2Si3O12(g)",1,1,"Ca3Al2Si3O12(g)",1,1,"Mg4Si4O12(g)",1,2,1);
    solution opxM1("MgSiO3(o)",0.5,1, "FeSiO3(o)",0.5,2, "MgAl2SiO6(o)",1,1,1);
    solution cpxSi("CaMgSi2O6(c)",1,1,"CaFeSi2O6(c)",1,1,"CaAl2SiO6(c)",1,2,1);
    
    ssref[8] = opxM2;     ssref[9] = cpxM1;  ssref[10] = garnetX;
    ssref[11] = garnetY;  ssref[12]= opxM1;  ssref[13] = cpxSi;

    /* add relevant solid solution to "ss_list" from "ssref" */
    for (int k=0; k<(int)ssref.size(); k++){
        for (int j=0; j<(int)n_in.size(); j++){
            string _s = list[j].getMoleculeName();
            if (ssref[k].is_insystem(_s)){
                ss_list.push_back(ssref[k]);
                continue;
            }
        }
    }
    
    /* create matrix */
    int n_ss = (int)ss_list.size(),  m = (int)n_in.size();
    ss_matrix.resize(0.0, m, n_ss);
    for (int j=0; j<m; j++){
        string _s = list[j].getMoleculeName();
        for (int k=0; k<n_ss; k++){
            if (ss_list[k].is_insystem(_s)){ ss_matrix[j][k] = 1; }
        }
    }
}
