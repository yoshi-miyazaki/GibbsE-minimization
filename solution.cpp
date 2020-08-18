/*
  solid_solution.cpp
  ... Gibbs free energy description of solid solution
  
  Aug 2, 2017 by Yoshinori
 */

#include "solution.h"

/*--------------------------------------------------------------------------
 // class: melt
 ---------------------------------------------------------------------------*/
melt_mixing::melt_mixing(string s1, string s2, string s3){
    /* store end members into the object */
    end_members.resize(3);
    nmole.resize(0.,3);   nmole = 0.0;
    Wmix.resize(0.,2);    Wmix  = 0.0;
    end_members[0] = s1;  end_members[1] = s2;  end_members[2] = s3;    
}
melt_mixing::melt_mixing(string s1, string s2, string s3, tensor1d<double> _v){
    /* store end members into the object */
    end_members.resize(3);
    nmole.resize(0.,3);   nmole = 0.0;
    end_members[0] = s1;  end_members[1] = s2;  end_members[2] = s3;
    
    /* Margules paramter */
    if (_v.size() == 2){ Wmix = _v; }
    else { cout << "wrong vector size - Wmix..." << endl; exit(1); }
}
void melt_mixing::add_endmember(string s1, string s2, string s3){
    /* store end members into the object */
    end_members.resize(3);
    nmole.resize(3);  nmole = 0.0;
    end_members[0] = s1;  end_members[1] = s2;  end_members[2] = s3;
}
void melt_mixing::add_Margules(double T, double P, tensor1d<double>& v){
    if (v.size()!= 6){ cout << "in add_Margules(): wrong v.size()... " << endl; exit(2);}
    
    /* calc Margules parameter */
    Wmix.resize(0.0,2);  Wmix = 0.0;
    Wmix[0] = v[0] + v[1]*P - v[2]*(T-1500);
    Wmix[1] = v[3] + v[4]*P - v[5]*(T-1500);
    Wmix /= (R*T);
}
bool melt_mixing::is_insystem(string _s){
    /* return whether _s is in this melt-mixing */
    bool is_in = 0;
    int no_endmember = end_members.size();
    for (int i=0; i<no_endmember; i++){
        if (_s == end_members[i]){ is_in = 1; break; }
    }
    return is_in;
}
void melt_mixing::add_mole(string _s, double n){
    /* store mole amount of end member _s */
    int no_endmember = end_members.size();
    for (int i=0; i<no_endmember; i++){
        if (_s == end_members[i]){ nmole[i] = n; break; }
    }
}
double melt_mixing::tot_mole(){
    /* return total amount in this solid-solution system */
    int    no_endmember = end_members.size();
    double Ntot = 0.0;
    for (int i=0; i<no_endmember; i++){ Ntot += nmole[i]; }
    return Ntot;
}
double melt_mixing::activity(string _s){
    /* return activity (RT*ln(gi)) of _s */    
    /* calculate total mole in this melt-mixing system */
    double Ntot = tot_mole();
    
    /* molar proportions in liquid phase */
    double xMg = 0., xFe = 0., xSi = 0.;
    for (int i=0; i<nmole.size(); i++){
        string _t = end_members[i];
        if     (_t == "MgO(l)"){ xMg = nmole[i]/Ntot; }
        else if(_t == "FeO(l)"){ xFe = nmole[i]/Ntot; }
        else if(_t =="SiO2(l)"){ xSi = nmole[i]/Ntot; }
    }
    
    /* calc activity */
    double act = 0.0;
    if      (_s == "MgO(l)"){ act = ((1-xMg)*Wmix[0] - xFe*Wmix[1])*xSi;  }
    else if (_s == "FeO(l)"){ act = ((1-xFe)*Wmix[1] - xMg*Wmix[0])*xSi;  }
    else if (_s =="SiO2(l)"){ act = (xMg*Wmix[0] + xFe*Wmix[1])*(xMg+xFe);}
    
    // cout << "act: " << act << endl;
    return act;
}
/*--------------------------------------------------------------------------
 // class: solution
 ---------------------------------------------------------------------------*/
solution::solution(string s1, string s2, int n){
    /* store end members (s1 & s2) into the object */
    end_members.resize(2);
    contribution.resize(1,2);  group.resize(0,2);
    nmole.resize(0.,2);
    end_members[0] = s1;  end_members[1] = s2;
    group[0]       = 1;   group[2]       = 2;
    
    no_sites = n;
}
solution::solution(string s1, string s2, string s3, int n){
    /* store end members (s1, s2 & s3) into the object */
    end_members.resize(3);
    contribution.resize(1,3);  group.resize(0,3);
    nmole.resize(0.,3);
    end_members[0] = s1;  end_members[1] = s2;  end_members[2] = s3;
    group[0]       = 1;   group[1]       = 2;   group[2]       = 3;
    
    no_sites = n;
}
solution::solution(string s1, double c1, int g1, string s2, double c2, int g2, string s3, double c3, int g3, int n){
    /* store end members (s1, s2 & s3) into the object */
    end_members.resize(3);
    contribution.resize(1,3);  group.resize(0,3);
    nmole.resize(0.,3);
    end_members[0]  = s1;  end_members[1]  = s2;  end_members[2]  = s3;
    contribution[0] = c1;  contribution[1] = c2;  contribution[2] = c3;
    group[0]        = g1;  group[1]        = g2;  group[2]        = g3;
    
    no_sites = n;
}
solution::solution(string s1, string s2, string s3, string s4, int n){
    /* store end members (s1, s2, s3 & s4) into the object */
    end_members.resize(4);
    contribution.resize(1,4);  group.resize(0,4);
    nmole.resize(0.,4);
    end_members[0] = s1;  end_members[1] = s2;  end_members[2] = s3;  end_members[3] = s4;
    group[0]       = 1;   group[1]       = 2;   group[2]       = 3;   group[3]       = 4;
    no_sites = n;
}
solution::solution(string s1, double c1, int g1, string s2, double c2, int g2, string s3, double c3, int g3, string s4, double c4, int g4, int n){
    /* store end members (s1, s2, s3 & s4) into the object */
    end_members.resize(4);
    contribution.resize(1,4);  group.resize(0,4);
    nmole.resize(0.,4);
    end_members[0]  = s1;  end_members[1]  = s2;  end_members[2]  = s3;  end_members[3]  = s4;
    contribution[0] = c1;  contribution[1] = c2;  contribution[2] = c3;  contribution[3] = c4;
    group[0]        = g1;  group[1]        = g2;  group[2]        = g3;  group[3]        = g4;
    
    no_sites = n;
}
bool solution::is_insystem(string _s){
    /* return whether _s is in this solid solution system or not... */
    bool is_in = 0;
    int no_endmember = end_members.size();
    for (int i=0; i<no_endmember; i++){
        if (_s == end_members[i]){ is_in = 1; break; }
    }
    return is_in;
}
void solution::add_mole(string _s, double n){
    /* store mole amount of end member _s */
    int no_endmember = end_members.size();
    for (int i=0; i<no_endmember; i++){
        if (_s == end_members[i]){ nmole[i] = n; break; }
    }
}
double solution::tot_mole(){
    /* return total amount in this solid-solution system */
    int    no_endmember = end_members.size();
    double Ntot = 0.0;
    for (int i=0; i<no_endmember; i++){ Ntot += contribution[i]*nmole[i]; }
    return Ntot;
}
double solution::config(string _s, double nthis){
    /* return configuration entropy 
       ... nthis = 0.0 => just output the log(xi) (activity) term
           nthis > 0   => specify the amount of [_s] 
                          even though it may have a diffrent amount 
                          (especially used when nthis is tiny (non0)  */
    int no_endmember = end_members.size();   int no_s;
    for (int i=0; i<no_endmember; i++){
        if (_s == end_members[i]){ no_s = i; break; }
    }
    // if (nmole[no_s] == 0){ return 0; }   /* if _s is not presenet, config S = 0 */
    
    /* calculate total mole in this solid-solution system */
    /* configuration entropy would be the log of concentration of _s in the cite.
       ... but for example in pyrope-majorite case, both have Mg for X cite
           in those case, group these two up, and divide it by the whole garnet amount. */
    double Ntot = tot_mole();
    double Ngroup = 0.0;
    for (int i=0; i<no_endmember; i++){
        // cout << "spec " << i << " - " << _s << " : " << no_s << " : " << Ngroup << endl;
        // cout << "group: " << group[i] << " - this: " << group[no_s] << endl;
        if (i == no_s && nthis > 0.0 )   { Ngroup += contribution[i]*nthis;    }
        else if (group[i] == group[no_s]){ Ngroup += contribution[i]*nmole[i]; }
    }    
    // cout << "config: " << log(Ngroup/Ntot) << " , " << nthis << " , group " << Ngroup << endl;
    return no_sites*contribution[no_s]*log(Ngroup/Ntot);
}
