/*
  solution.h
  ... part of CGgibbsmin.cpp
  
  created by Yoshi Miyazaki on 2018/07/04
 */

#ifndef GIBBS_solution_h
#define GIBBS_solution_h

#include <iostream>
#include "./tensor.h"
#include "./matrix.h"
#include "./const.h"

/*--------------------------------------------------------------------------
  define `melt_mixing` class
 ---------------------------------------------------------------------------*/
class melt_mixing{
 public:
    melt_mixing(string s1 = "", string s2 = "", string s3 = "");
    melt_mixing(string s1, string s2, string s3, tensor1d<double>);
    
    void    add_endmember(string, string, string);
    void    add_Margules(double, double, tensor1d<double>&);
    
    bool    is_insystem(string);
    void    reset(){ nmole = 0; };
    void    add_mole(string, double);   /* store mole amount of end member */
    double  tot_mole();                 /* return total amount in this solid-solution */
    double  activity(string);
 private:
    tensor1d<string>    end_members;
    tensor1d<double>    Wmix;
    tensor1d<double>    nmole;
};

/*--------------------------------------------------------------------------
  define `solution` class
 ---------------------------------------------------------------------------*/
class solution{
 public:
    solution(string s1 = "", string s2 = "", int no_sites=0);
    solution(string s1, string s2, string s3, int no_sites);
    solution(string s1, string s2, string s3, string s4, int no_sites);
    solution(string, double, int, string, double, int, string, double, int, int);
    solution(string, double, int, string, double, int, string, double, int, string, double, int, int);
    
    int     get_nosites(){ return no_sites; }
    bool    is_insystem(string);
    void    reset(){ nmole = 0; };
    void    add_mole(string, double);   /* store mole amount of end member */
    double  tot_mole();                 /* return total amount in this solid-solution */
    double  config(string, double);     /* return configuration entropy               */
 private:
    tensor1d<string>    end_members;
    int                 no_sites;
     
    tensor1d<int>       group;
    tensor1d<double>    contribution;
    tensor1d<double>    nmole;
};

#endif
