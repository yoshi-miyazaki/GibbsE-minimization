/*
  gibbsminCG.cpp
  Gibbs energy minimization by conjugate gradient method
  
  by Yoshi Miyazaki on August 4, 2015
  All rights reserved.
*/

/*--------------------------------------------------------------------------
 // Def of class: gibbsminCG
 ---> Refer to ./CGgibbsmin.h for explanation
 ---------------------------------------------------------------------------*/

#include "CGgibbsmin.h"

#include <float.h>

gibbsminCG::gibbsminCG(MoleculeData_G& moleculelist, tensor1d<double>& n_starting){
    tensor1d<double>  W0(0.0, 6);
    // W0[0]=-125000;  W0[1]=0.088e-06;  W0[2]=30;
    // W0[3]= -45200;  W0[4]=0.190e-06;  W0[5]=0.;
    W0[0]=-126705;  W0[1]=1.353e-06;  W0[2]=3.44128;
    W0[3]=-99917.;  W0[4]=6.329e-07;  W0[5]=52.9655;
    
    *this = gibbsminCG(moleculelist, n_starting, W0);
}
gibbsminCG::gibbsminCG(MoleculeData_G& moleculelist, tensor1d<double>& n_starting,
                       tensor1d<double> W0){
    /* store basic information */
    list = moleculelist;   //list.rm_nonexisting(n_starting);
    
    /* pressure and temperature + thermodynamic data  */
    T = list[0].getT();   P = list[0].getP();
    gibbse = getGibbslist(list);      /* Gibbs energy (converted to c_j)        */
    phase  = getphaselist(list);      /* phase                                  */
    
    /* non-ideal melt and solid-solution */
    meltSi.add_endmember("MgO(l)","FeO(l)","SiO2(l)");
    meltSi.add_Margules(T,P,W0);
    create_sslist();
    
    /* mass-balance matrix */
    massBalance  massbo(list);
    massm  = massbo.get_massBalance();
    atomic = massbo.get_atomic();
    
    /* conjugate gradient */
    n_in = n_starting;                      /* starting vector for conjugate gradient */
    conjugateG();
}
void gibbsminCG::conjugateG(){
    int m = n_in.size();    /* size of system */
    // srand((unsigned int)time(NULL));
    
    /*cout << "mass balance" << endl;
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
    int  count = 0, count_cg = 0, bound = 10, CMAX = 1000;
    double eps = 1e-6;
    
    while (count < CMAX){
        count++;  count_cg++;
        tensor1d<double> n_save = n_old;
	
        /* store results from previous step. */
        n_old = n;  dn_old = dn;  oldG = G;  exist_old = exist;
	
        /* delete un-resolvable numbers */
        double n_unres = n_old.max()*(LIM*0.1);
        for (int i=0; i<m; i++){ if(n_old[i] < n_unres && n_old[i] > 0){ n_old[i] = 0.0;}}
        for (int i=0; i<m; i++){ if (std::isnan(n_old[i])){ cout << "error, nan in : " << i << " at step " << count <<endl; }}
	
        /* calculate gradient */
        exist = 1;    dn = grad(n_old,exist)*(-1.0);
        /* trick to converge faster. if gas/liq mole i = 0 => exist[i] = 0 */
        if (count <= bound){ comp = 1; }
        if (count > bound && comp == 0){ zerogas(n_old,dn,exist); }
        avoidneg(n_old,dn,exist);  /* project to balance mass & delete unnessary spec */
	
        /* conjugate gradient
           ... when the dimension changeed, conguate => steepest descent */
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
        if (cg2small || cg.norm() < eps){           /* if cg is too small (cg2small) */
            if (comp == 0){ comp = 1; continue; }}  /* AND comp = 0, reset the search and
                                                       calculate cg in complete dimension */
	
        /* termination criteria => |dn| < eps for M consecutive times */
        if (dn.norm() < eps && comp > 0){/* but only if dn is calced in complete dimension
                                            (which means comp == 1). */
            int M = m*10;
            for (comp=1; comp<M; comp++){   /* when i = i, comp = (i+1) */
                exist = 1;
		dn = grad(n_old,exist)*(-1.0);

		if (comp==1){
		    dn = grad(n_old + n_in*DBL_MIN, exist)*(-1.0);
		}

		/*
		for (int j=0; j<m; j++){
		    cout << setw(15) << list[j].getMoleculeName() << " - " << setw(10) << setprecision(10) << n_old[j] << "  -  " << setw(20) << dn[j] << " - " << exist[j] << endl;
		}
		*/
		
                avoidneg(n_old,dn,exist);

		/*
		for (int j=0; j<m; j++){
		    cout << setw(15) << list[j].getMoleculeName() << " - " << setw(10) << setprecision(10) << n_old[j] << "  -  " << setw(20) << dn[j] << " - " << exist[j] << endl;
		}
		*/
		
                /* check whether any direction is heading towards smaller G */
                double  dc = 1;   n = shrink(n_old,dn,dc);
                bool small = ck_smaller(n_old,n,dn,exist);
		
		// cout << dn.norm() << "\t" << small << endl;
                // if (dn.norm() > eps && small){ cg = dn; break; }
		if (dn.norm() > eps){ cg = dn; break; }
            }
	    if (comp == M){ cg = dn;  break; }
        }
        comp = 0; /* reset comp here. */
	
        /* check the sign of cg. 
           ... even if n[i] = 0 and dn[i] >=0, cg[i] < 0 may occur.
           in such cases, reset beta = 0 => cg = dn */
        for (int i=0; i<m; i++){
            if (n_old[i]==0.0 && cg[i]< 0){ beta = 0.0;  cg = dn; break; }	
        }
	
        /* next point / modify(*dn) : factor to adjust neg mole   */
        double mod;
        n = shrink(n_old,cg,mod);              /* shrink arrow to avoid neg mole       */
        tensor1d<double> nright = n;
        bisection(n_old,n,cg,mod,exist,sd,up); /* bisection or golden section search   */
        nbest  = n;  G = gibbsE(nbest);        /* when dG = 0 is not found, sd will be
                                                  updated to 1, and the beta next step 
                                                  is set to 0. (sd = steepest descent) */
	
        if (count > 30000){
            tensor1d<double> gc = grad(nbest,exist);
            cout << endl;
            cout << count << " ..c = " << comp << ", up: " << up << ", sd?: " << sd << endl;
            cout << "G - oldG : " << fixed << G-oldG << ", beta: " << beta <<  endl;
            for (int i=0; i<m; i++){
                cout << setw(15) << list[i].getMoleculeName();
                cout << " -old: " << setw(10) << setprecision(7) << scientific << n_old[i];
                cout << " -new: " << n[i]  << "  - shrinked: " << nright[i];
                cout << " - dn = " << setw(10) << fixed << dn[i] << " - cg: " << cg[i];
                cout << " ex= " <<  fixed << exist[i] << " - mu = " << gc[i] << "\t" << gibbse[i] << endl;
            }
            cout << "done. " << cg.norm() << endl << endl;
        }
    }
    
    
    /* output the final result. 
    if (count > 2990 or 0){
        cout << "T = " << T << " [K]" << " , |dG/dn| = " << cg*cg;
        cout << " ,c = " << count << " ,ex= " << exist.sum() << " ,comp= " << comp << endl;
	}
    dn = grad(nbest,exist);  tensor1d<double> gc = grad(nbest,exist);
    avoidneg(n_old,dn,exist); 
    
    cout << "T = "<< T << " [K]" << " , |dG/dn| = "<< dn*dn<<" , G = "<< G<<endl;
    cout << "-               n:          cg:              dn:            exist:" << endl;
    for (int j=0; j<m; j++){
        cout << setw(17) << nbest[j] << "  " << setw(16) << cg[j] << "  " << setw(16) << gc[j] << " " << exist[j] << " - " << list[j].getMoleculeName() << endl;}
    cout << endl;  
    cout << "... " << P/1e9 << " GPa, c= " << count << endl;      */
    /*cout << T << "\t" << P/1e9 << "\t" << G << "\t";
      for (int j=0; j<(int)nbest.size(); j++){ cout << nbest[j] << "\t"; } 
      cout << count << endl; */

    /* elimenate very small amount */
    for (int i=0; i<m; i++){ if (n[i] < 1e-300){ n[i] = 0.; } }

    dn_last = dn;
    return;
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
    int m = n_in.size();            /* size of system               */
    int n_mixing = ss_list.size();  /* no. of solid solution models */
    
    /* delete all mole & configuration entropy in each "solid_solution" obejct */
    meltSi.reset();
    for (int k=0; k<n_mixing; k++){ ss_list[k].reset(); }
    
    /* remove unresolvable number */
    for (int i=0; i<n.size(); i++){ if (n[i] > 0 && n[i] < DBL_MIN){ n[i] = 0.0; }}
    
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
                if (ss_list[k].is_insystem(_s)){ ss_list[k].add_mole(_s,n[j]); }
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
            } else {
                string _s = list[j].getMoleculeName();
                for (int k=0; k<n_mixing; k++){
                    if (ss_list[k].is_insystem(_s)){ dG+= n[j]* ss_list[k].config(_s, 0.0);} 
                }/* end of solid solution model search */
            }/* end of gas/liquid/ss search */
        }
        G += dG;  /* add to total Gibbs free energy */
    }
        
    return G;
}
tensor1d<double> gibbsminCG::grad(tensor1d<double> n, tensor1d<int> exist){
    int m = n_in.size();            /* size of system               */
    int n_mixing = ss_list.size();  /* no. of solid solution models */
    
    /* remove unresolvable number */
    for (int i=0; i<n.size(); i++){ if (n[i] > 0 && n[i] < DBL_MIN){ n[i] = 0.0; }}
    
    /* delete all mole & configuration entropy in each "solid_solution" obejct */
    meltSi.reset();
    for (int k=0; k<n_mixing; k++){ ss_list[k].reset(); }
    
    /* calculate mixing of gas, liquid, and solid-solution (ss) */
    double Ngas = 0.0, Nliq = 0.0;  int type_liq = 0;
    tensor1d<int>    type_ss(0,n_mixing); /* store number of species in each ss system */ 
    tensor1d<double> Nss(0.0,  n_mixing); /* store total mole of each ss system        */
    for (int j=0; j<m; j++){
	if      (exist[j] == 0){ continue; }
        if      (phase[j] == 0){ Ngas += n[j]; }
        else if (phase[j] == 1){ Nliq += n[j];  type_liq += 1;
            string _s = list[j].getMoleculeName();
            meltSi.add_mole(_s, n[j]);
        } else {
            string _s = list[j].getMoleculeName();
            for (int k=0; k<n_mixing; k++){
                if (ss_list[k].is_insystem(_s)){
		    type_ss[k] += 1;
		    ss_list[k].add_mole(_s,n[j]);
		    break;
		}
	    }
        }
    }
    
    /* for gas, liquid, and solid-solution,
       when their total amount is 0, no lowering through activity happens. 
       A different phase may be showing a lower chemical potential, but
       these phases may be more stable when mixing is included.
       
       To avoid this, the code assigns a very small amount to calculate the gradient. */
    
    if (Nliq == 0){ // for liquid phase
        Nliq = 1.0e-15;  double Nrm = Nliq;  int type_rm = 0;
        for (int j=0; j<m; j++){
            if (phase[j] == 1 && exist[j] == 1){
                type_rm++;
                if (type_rm != type_liq){ n[j] = Nrm*(double)rand()/RAND_MAX; }
                else                    { n[j] = Nrm; }
                Nrm -= n[j];
		
		// update moles in the melt-activity model
                string _s = list[j].getMoleculeName();
                meltSi.add_mole(_s, n[j]);
            }
        }
    }
    
    for (int k=0; k<n_mixing; k++){ // for solid solutions
	if (k==0)
	if (ss_list[k].tot_mole() == 0){
	    double Ns = 1.0e-15; double Nrm = Ns; int type_rm = 0;
	    for (int j=0; j<m; j++){
		string _s = list[j].getMoleculeName();
		if (ss_list[k].is_insystem(_s) && exist[j] == 1){
		    type_rm++;
		    if (type_rm != type_ss[k]){
			// use logit distribution to avoid being stuck at n = 0
			double rn = -8. + 16.*(double)rand()/RAND_MAX;
			n[j] = Nrm/(1.+exp(-rn));
			// if (k==0){ cout << "* " << scientific << rn << endl; }
		    }
		    else { n[j] = Nrm; }
		    Nrm -= n[j];
		    
		    // update moles in the solid solution model
		    ss_list[k].add_mole(_s, n[j]);
		}
	    }
	    // if (k==0 and exist.sum() == 7){ cout << k << "\t" << n[2] << "\t" << n[3] << "\t" << n[2]/n[3] << endl; }
	}
    }

    
    /* obtain Mg, Fe, Si element mole   
    Matrix<double>    trans = massm.transpose();
    tensor1d<double>  q_elm = trans*n;
    
    int eMg = -1, eSi = -1, eFe = -1;
    for (int l=0; l<q_elm.size(); l++){
        if (atomic[l] == 12){ eMg = l; }
        if (atomic[l] == 14){ eSi = l; }
        if (atomic[l] == 26){ eFe = l; }
    }
    */
    
    /* calculate dG gradient in a space*/
    tensor1d<double> dG(m);
    double           non0 = abs(n.max()) * LIM;
    for (int j=0; j<m; j++){
        string _s = list[j].getMoleculeName();
        dG[j] = gibbse[j];
	if (exist[j] == 0) continue;
	
        /* for gas, liquid, and solid-solution, assume ideal mixing
           ... chemical potential should be modified as mu = mu0 + RT*log(concentration)
           BUT for very small concentration, smaller than the numerical digit 
           than C++ could handle, assume non0 mole instead of actual amount.
           This would prevent large gradient value and make converge faster */
        if (phase[j] == 0){
            if (n[j] > 0)    { dG[j] = gibbse[j] + log(n[j] / Ngas);}
            else if(Ngas > 0){ dG[j] = gibbse[j] + log(non0 / Ngas);}
        }else if(phase[j] == 1){
            double dGmix = meltSi.activity(_s);
            if (n[j] > 0)    { dG[j] = gibbse[j] + log(n[j] / Nliq) + dGmix;}
	    else if(Nliq > 0){ dG[j] = gibbse[j] + log(non0 / Nliq) + dGmix;}
        }else {
            dG[j] = gibbse[j];
            for (int k=0; k<n_mixing; k++){
                if (ss_list[k].is_insystem(_s)){
                    double Ntot  = ss_list[k].tot_mole();
                    // cout << " -- " << _s << " , " << ss_list[k].config(_s, 0.0) << endl;
                    if  (n[j] > 0){ dG[j] += ss_list[k].config(_s, 0.0); }
                    else if(Ntot > 0){ dG[j] += ss_list[k].config(_s,non0); }
                    // cout<<"n0: "<< _s << " = " << scientific <<  ss_list[k].config(_s,non0) << " , non0 = " << non0 << " , nmax = " << abs(n.max()) << endl; }
                    break;
                }
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
    
    /* if (isnan(dn[0])){
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
    int m = n_in.size();    /* size of system */
    int n_mixing = ss_list.size();
    /* delete all mole & configuration entropy in each "solid_solution" obejct */
    for (int k=0; k<n_mixing; k++){ ss_list[k].reset(); }
    
    /* calculate mixing of gas, liquid, and solid-solution (ss) */
    double Ngas = 0.0, Nliq = 0.0;
    tensor1d<double> Nss(0.0, n_mixing); /* store total mole of each ss system */
    for (int j=0; j<m; j++){
        if      (phase[j] == 0){ Ngas += n_old[j]; }
        else if (phase[j] == 1){ Nliq += n_old[j]; }
        else { string _s = list[j].getMoleculeName();
            for (int k=0; k<n_mixing; k++){
                if (ss_list[k].is_insystem(_s)){ ss_list[k].add_mole(_s,n_old[j]); }}
        }
    }
    
    for (int j=0; j<m; j++){
        if(n_old[j] == 0){
            if      (phase[j]==0){ dn[j]=0.0; exist[j]=0; }
            else if (phase[j]==1){ dn[j]=0.0; exist[j]=0; }
            else { string _s = list[j].getMoleculeName();
                for (int k=0; k<n_mixing; k++){
                    if (ss_list[k].is_insystem(_s)){
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
    /* 
       modify `dn` to avoid searching 
     */
    int m = n_in.size();    /* size of system */
        
    tensor1d<double> dt = dn;
    int count = 0, ifneg = -2;
    
    while (ifneg < 0){
	count++;
	dt = grad(n_old,exist) * (-1.0);

	for (int j=0; j<m; j++){ if(exist[j] == 0){dt[j] = 0.0;} }
	mass_balance(dt,exist);         /* project to meet mass balance eq.  */	

	/*
	if (exist.sum() == 7 and n_old[2] == 0){
	    for (int j=0; j<m; j++){
		cout << setw(15) << list[j].getMoleculeName() << " + " << setw(10) << setprecision(10) << fixed << n_old[j] << "  -  " << setw(20) << dt[j] << " - " << exist[j] << endl;
	    }
	    cout << endl;
	}
	*/
	
	/* check whether it satisfies dt[j] >= 0 when n_old[j] = 0.0 */
	ifneg = 1;
	int s = rand()%m;
	for (int j=0; j<m; j++){
	    if (n_old[s] == 0 && dt[s] < 0){ ifneg = -1; exist[s] = 0; break;}
	    s++; s=s%m;
	}
    }

    /*
    for (int j=0; j<m; j++){
	cout << setw(15) << list[j].getMoleculeName() << " + " << setw(10) << setprecision(10) << fixed << n_old[j] << "  -  " << setw(20) << dt[j] << " - " << exist[j] << endl;
    }
    cout << endl;
    */
    dn = dt;

}
tensor1d<double> gibbsminCG::shrink(tensor1d<double>& n_old, tensor1d<double>& dn, double& mod){
    /* adjust dn so that n = n_old+dn have non-negative amount */
    int m = n_in.size();    /* size of system */
    tensor1d<double> n  = n_old +dn;
    tensor1d<double> dt = dn;
    
    /* shrink arrow */
    double temp;
    mod = 1.0;   dore = -1;     
    for (int j=0; j<m; j++){
	if (n_old[j] > 0 && n[j] < 0){
	    temp = abs(n_old[j])/(abs(n_old[j])-n[j]);
	    if (temp < mod){mod = temp; dore = j;}
	}
    }
    dt *= mod;    n = n_old +dt;
    if (dore != -1){ n[dore] = 0.0;}    /* assure n[dore] = 0 */
    
    /* reassure that no negative mole exists */
    double negLIM = n_old.max()*(-1.0*LIM);
    for (int j=0; j<m; j++){
	if (negLIM < n[j] && n[j] < 0){ n[j] = 0.0; }
	// if (n[j] < 0){ cout << "negative mole >_< : " << j << " is " << scientific << n[j] << endl;}
    }
    
    return n;
}
double gibbsminCG::dG_direction(tensor1d<double>& n, tensor1d<int>& exist, tensor1d<double> cg){
    /* calculate change in GFE along cg direction. dG/da = (nabla G)*cg */
    tensor1d<double> n_pert = perturb_zero(n,cg);
    tensor1d<double> nabG   = grad(n_pert,exist); // mass_balance(nabG,exist);
    //for (int i=0; i<(int)nabG.size(); i++){cout << n[i] << " - " << nabG[i] << endl; }
    
    double dGda = nabG*cg;
    if (std::isnan(dGda)){
        cout << "nan... @" << T << "K, P: " << P/1e9 << endl;
        for (int j=0; j<n.size(); j++){
            cout << j <<" - "<< n_pert[j] <<" - cg: "<< cg[j] <<" , grad: "<< nabG[j] << endl;
        }
    }
    return dGda;
}
bool gibbsminCG::ck_smaller(tensor1d<double>& n_old, tensor1d<double>& n, tensor1d<double>& dn, tensor1d<int>& exist){
    /* judge whether direction dn is likely to lower GFE 
       ... check using the change in sign of r = (nabla G)*(dn) */
    bool   Gsmaller = 0;    
    double oldG = gibbsE(n_old), G = gibbsE(n);
    if (oldG > G){ Gsmaller = 1; }
    else { /* if bisection search is possible, 1 */
        double oldr = dG_direction(n_old,exist,dn);
        double    r = dG_direction(n,    exist,dn*(-1.))*(-1.);
        if (oldr*r < 0){ Gsmaller = 1;}
        //cout << "r = " << r << " , oldr = " << oldr << " Gsm? = " << Gsmaller << endl;
	
        // double refr = dG_direction(n_old,exist,dn);
    }
    
    //cout << "-- G = " << G << " , oldG = " << oldG << endl;
    //for (int i=0; i<(int)n.size(); i++){ cout << n_old[i] << " - " << n[i] << endl; }
    
    return Gsmaller;
}
tensor1d<double> gibbsminCG::perturb_zero(tensor1d<double>& n, tensor1d<double>& dn){
    tensor1d<double> n_pert = n;
    
    double mindn = 0.0; /* perturb zero-mole species to non0 if dn(i) != 0 */
    for (int i=0; i<n.size(); i++){ if (n[i] == 0 && dn[i] > 0){ mindn = maxv(mindn,dn[i]); }}
    if (mindn == 0.0){ return n_pert; }
    
    /* reject any negative mole */
    double  non0 = DBL_MIN*1e5;
    n_pert += dn*(non0/mindn);
    // cout << "ratio: " << scientific << non0/mindn << endl;

    double mod = 1.0;
    for (int i=0; i<n.size(); i++){
	if (n_pert[i] < 0){ mod = minv(mod, n[i]/(n[i]-n_pert[i]));}
	// cout << "pert: " << (non0/mindn)*mod << " , n-n_pert : " << n[i]-n_pert[i] << endl;
    }
    n_pert = n + dn*(non0/mindn)*mod*(1-DBL_MIN);
    
    return n_pert;
}
void gibbsminCG::bisection(tensor1d<double> n_old, tensor1d<double>& n, tensor1d<double>& cg, double mod, tensor1d<int>& exist, int& sd, int& up){ /* cg: direction (before shrink) */
    int m = n_in.size();     /* size of system  */
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
    tensor1d<double> g = grad(n,exist);   mass_balance(g,exist);
    double  oldr = dG_direction(n_old,exist, cg);
    double     r = dG_direction(n,    exist, cg*(-1.))*(-1.);
    
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
            if (toosmall == 1){ // cout << "too small: elements:" << endl;
                /*for (int k=0; k<e; k++){
                  cout << "    --- compare: " << scientific << abs(q_old[k]-q[k]) << " - ";
                  cout <<  q_old[k]*LIM << endl;
                  }*/
                break;
            }
	    
            /* break loop if any of the species are negative */
            int ifneg = 1;
            for (int j=0; j<m; j++){ if(n[j] < 0){ ifneg = -1; }}
            if (ifneg < 0){
                n = shrink(n_old, cgs, mod);
                r = dG_direction(n,exist,cg*(-1.))*(-1.);   break;
            }
	    
            /* for the next loop */
            r = dG_direction(n,exist,cg*(-1.))*(-1.);  /* renew the value of r */
        }
    }
    // cout << " -- pert, r signs: left: " << oldr << " , right: " << r << endl;
    /* (just in case:) check if (n) diverged by doubling (cg) */
    for (int j=0; j<m; j++){
        if (std::isnan(n[j])){ cout << list[j].getMoleculeName() << " has diverged during cg amplification... " << n[j] << endl; exit(10);}}
    
    /* step-3. actual LINE SEARCH : golden seciton search or biseciton
       ... res = absnon0 min value of vector |n-n_old| */
    tensor1d<double> nleft = n_old, nright = n, ds;
    double res  = 2.0,  oldres = 0.0;
    double oldG = gibbsE(n_old), G = gibbsE(n);
    
    bool yeah = 0;
    if (G - oldG > 128 and G - oldG < 130){ yeah = 1; }
    if (oldr*r > 0){ /* the same sign => golden section search */
        tensor1d<double> nA, nB, g_old, g;
        double GA, GB, tau = 0.61803398875;
        int refc = 0;
        while (oldres != res && refc < 40){
            refc++;
            ds = n - n_old;
            nA = n_old + ds*(1-tau);   GA = gibbsE(nA);
            nB = n_old + ds* tau;      GB = gibbsE(nB);
	    
            if (minv(GA,GB) < minv(oldG,G)){      /* switch */
                if (GA > GB){n_old= nA; oldG= GA; } /* next search b/w (nA   ,n ) */
                else        {n    = nB;    G= GB; } /*             b/w (n_old,nB) */
            } else if (minv(GA,GB) < G and yeah) {
		if (GA > GB){n    = nB;    G= GB; }
		else        {n    = nA;    G= GA; }
		}else{
		n_out = n;
		break;
	    }
	    
            /* calculate gradient */
            oldr = dG_direction(n_old,exist,cg);
            r    = dG_direction(n,    exist,cg);
            if(oldr*r < 0){ nright = n; break; } // cout << "gss: r- " << r << endl; break; }
	    
            /* termination condition*/
            ds = n - n_old;  n_out = (n_old+n)/2.0;
            oldres = res;    res = ds.absnon0min();
        }
        if (oldr*r > 0){ 	/* if minimization fails : next step => sd */
            if(G <= oldG){n_out= n;     sd= 1; up= 0;}
            else{
                if(up==0){n_out= n_old; sd= 1; up= 1;} /* if prev step: cg-> redo w/ sd */
                else     {n_out= n;     sd= 1; up= 0;} /*            sd-> accept G>oldG */
            }
        }
    }

    if (oldr*r < 0){ /* different signs => bisection search */
        tensor1d<double> nA, gA;
        double rl = oldr, rA, rr = r;
	
        int crep = 0, cmax = 70;
        while (oldres != res and crep < cmax){
            /* refc > 40 should be unnecessary... */
            crep++;
            ds = n - n_old;   nA = n_old + ds*0.5;
            rA = dG_direction(nA,exist,cg);
            // cout << "rA = " << rA << endl;
	    
            if (rA == 0)  {n = nA;  n_out = n;  break;}
            if (rA*rl < 0){n = nA;  n_out = n_old;  rr = rA;}
            else      {n_old = nA;  n_out = n;      rl = rA;}
	    
            /* termination condition*/
            ds  = n-n_old;
            oldres = res;  res = ds.absnon0min();  /* termination condition */
        }
        up = 0;

	// cout << "bisec, old= " << scientific << rl << " , r= " << rr << " , G-old= " << G-oldG << "\t" << fixed << crep << endl;

	nA = n_old + ds*0.5;
	rA = dG_direction(nA,exist,cg);
	ds = nA - n_old;	
	
        /* if rA is not even close to 0 (= when bisection failes) */
        if (abs(rA)>1){
            // cout << " bisection fail oldr = " << oldr << " , rA = " << rA << " , r = " << r << endl;
            // oldr & rA same sign => real solution is r side
            if (oldr*rA > 0){ n_out = n; }
            else            { n_out = n_old; }
        }

        /*cout << "match, refc:" << refc << ", oldres: " << scientific << oldres << " res = " << res << " , rA = " << rA << endl;
          tensor1d<double> mcg = cg*(-1.0);
          tensor1d<double> n_pl = perturb_zero(nleft,cg), n_pr = perturb_zero(nright,mcg);
          double r1 = dG_direction(n_pl,exist,cg), r2 = dG_direction(nright,exist,cg*(-1.))*(-1.);
          double r3 = dG_direction(nA,exist,cg);
          cout << "rl = " << r1 << " , rA = " << rA << " , rr = " << r2 << endl;
          if (nleft == n_old || nright == n){
          sd = 1;
          for (int i=0; i<ds.size(); i++){
          cout << " - " << setw(12) << list[i].getMoleculeName() << " = " << ds[i] << " , lef = " << nleft[i] << " , lp = " << n_pl[i] << " , A = " << nA[i] << " , rp = " << n_pr[i] << " , right = " << nright[i] << endl;
          }
          for (int i=0; i<ds.size(); i++){
          cout << " - " << setw(12) << list[i].getMoleculeName() << " = " << setw(15) << cg[i] << " , lef = " << n_old[i] <<  " , out = " << n_out[i] << " , rig = " << n[i] << endl;
          }
          }*/
    }
    
    n = n_out;
    if(std::isnan(n[0])){ cout << "NAN happened in bisection." << endl; exit(17); }
}
/* output result */
void gibbsminCG::result(double T, double P, string& filename){ /* T in K, P in Pa */
    int m = n_in.size();     /* size of system */
    
    G = gibbsE(nbest);
    double mf = melt_vfrac();
    ofstream fout(filename, ios::app);
    fout << T << " " << P/1e9 << " " <<setprecision(10) << G;
    for (int j=0; j<m; j++){
        fout << " " << setprecision(10) << nbest[j];
    }
    fout << " " << mf << endl;
    fout.close();
    
    tensor1d<double> mu = grad(nbest,exist_saved);
    cout << "------------( Result )------------" << endl << endl;
    cout << "At T = " << T << " K. or " << T - 273.15 << " C." << endl;
    cout << "At P = " << P/1e9 << " GPa." << endl;
    cout << "Minimum Gibbs free energy is " << endl;
    cout << "  G = " << fixed << setprecision(12) << G << endl;
    cout << "Composition: " << endl;
    for (int j=0; j<m; j++){
        cout << setw(3) << j << setw(20) << list[j].getMoleculeName() << " phase: " << phase[j] << ",mol = " << fixed << setprecision(10) << nbest[j] << " , mu: = " << mu[j] << endl;
    }
    cout << endl;
}
double gibbsminCG::melt_vfrac(){
    int m = n_in.size();
    tensor1d<double> sol(0.0,m), liq(0.0,m), gas(0.0,m);
    for (int j=0; j<m; j++){
        if (phase[j]==2){ sol[j] = nbest[j];}
        if (phase[j]==1){ liq[j] = nbest[j];}
        if (phase[j]==0){ gas[j] = nbest[j];}
    }
    
    tensor1d<double> e_mel=(massm.transpose())*liq, e_sol=(massm.transpose())*sol;
    double v_mel = e_mel[0]*11.3 + e_mel[1]*22.67 + e_mel[2]*12.5;
    double v_sol =(e_sol[0]*11.3 + e_sol[1]*22.67 + e_sol[2]*12.5)*1.03;
    
    return v_mel/(v_mel+v_sol);
}
double gibbsminCG::melt_mfrac(){
    int m = n_in.size();
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
double gibbsminCG::meltfrac(double T, double P){
    double mfrac = melt_vfrac();
    return mfrac;
}
double gibbsminCG::meltfrac(double T, double P, string& filename){
    double mfrac = melt_vfrac();
    
    ofstream fout(filename, ios::out | ios::app);
    fout << T << "\t" << P << "\t" << mfrac <<  endl;
    fout.close();

    return mfrac;
}



/*--------------------------------------------------------------------------
 // Def of class: solid solution
 ---------------------------------------------------------------------------*/
void gibbsminCG::create_sslist(){
    ss_list.resize(5);
    
    /* define the solid-solution system. 
     ... when defining the object, specify <string> of end members, and constitution no.
      e.g.) for olivine, Mg2SiO4(f) will be 2, because of 2 Mg-sites. */
        
    solution     perovskite("MgSiO3(p)" , "FeSiO3(p)" , 1);
    ss_list[0] = perovskite;
    
    // solution     postperovskite("MgSiO3(v)" , "FeSiO3(v)" , 1);
    // ss_list[2] = postperovskite;
    
    solution     mgwustite("MgO(p)"   , "FeO(p)"    , 1);
    ss_list[1] = mgwustite;
    
    solution     olivine("Mg2SiO4(o)", "Fe2SiO4(o)", 2);
    ss_list[2] = olivine;

    solution     ringwoodite("Mg2SiO4(r)", "Fe2SiO4(r)", 2);
    ss_list[3] = ringwoodite;

    solution     opx("MgSiO3(o)" , "FeSiO3(o)" , 1);
    ss_list[4] = opx;
    // cout << "D " << ss_list.size() << " " << ss_list.capacity() << " " << sizeof(olivine) << flush;
    // cout << " " << sizeof(ss_list[0]) << " " << &ss_list[0] << " " << &ss_list << endl;
    
}
