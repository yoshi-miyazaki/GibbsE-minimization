//
//  matrix.h
//  Define Class for Vector & Matrix
//
//  Created by Yoshi Miyazaki on 2014/10/13.
//

#ifndef tensor_h
#define tensor_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <string>
#include <typeinfo>
#include <cstdlib>
#include <cmath>
using namespace std;

/* return square, max, min */
template<class T>
inline const T square(const T a){
    return a*a;
}
template<class T>
inline const T cube(const T a){
    return a*a*a;
}
template<class T>
inline const T cuatro(const T a){
    return a*a*a*a;
}
template<class T>
inline const T maxv(const T a, const T b){
  if (a > b) {return a;}
  else {return b;}
}
template<class T>
inline const T minv(const T a, const T b){
  if (a < b) {return a;}
  else {return b;}
}

/*----------------------------------------
 1d-vector class
 ---------------------------------------*/
template<class T>
class tensor1d : public vector<T>{
 private:
    int nrow;
 public:
    explicit tensor1d(size_t n=0);
    explicit tensor1d(const T&, size_t);    /* initialize to constant a */
    explicit tensor1d(const T*, size_t);    /*            to an array   */
    
    /* operator */
    tensor1d& operator=(const tensor1d&);
    tensor1d& operator=(const T&);
    const  bool operator==(const tensor1d& rhs) const;
    const  bool operator!=(const tensor1d& rhs) const{ return !(*this == rhs);};
    inline T&       operator[](const int i)       { return vector<T>::operator[](i); };
    inline const T& operator[](const int i) const { return vector<T>::operator[](i); };
    
    /* mathematical operation */
    const T    norm()    const;
    const T    sum()     const;
    const T    max()     const;
    const T    min()     const;
    const int  argmax()  const;
    const int  argmin()  const;
    const bool isnan()   const;
    const T    absnon0min()  const;
    
    tensor1d<T>  operator+ (const tensor1d<T>&) const;
    tensor1d<T>  operator- (const tensor1d<T>&) const;
    T            operator* (const tensor1d<T>&) const;
    //tensor1d<T>  operator* (const tensor1d<T>&) const;
    tensor1d<T>  operator/ (const tensor1d<T>&) const;
    tensor1d<T>  operator* (const T&) const;
    tensor1d<T>  operator/ (const T&) const;
    tensor1d<T>& operator+=(const tensor1d<T>&);
    tensor1d<T>& operator-=(const tensor1d<T>&);
    tensor1d<T>& operator*=(const T&);
    tensor1d<T>& operator/=(const T&);
    
    /* inherited from vector class */
    size_t   size() const           { return vector<T>::size();      };
    void     resize(size_t n)       { vector<T>::resize(n);          };
    void     resize(T var, size_t n){ vector<T>::resize(n, var);     };
    void     erase(int ir)          { vector<T>::erase(begin()+ir);  };
    
    typename vector<T>::const_iterator cbegin() const { return vector<T>::cbegin(); };
    typename vector<T>::const_iterator cend()   const { return vector<T>::cend(); };
    typename vector<T>::iterator begin() { return vector<T>::begin(); };
    typename vector<T>::iterator end()   { return vector<T>::end(); };
};

#include "tensor_s.h"

#endif
