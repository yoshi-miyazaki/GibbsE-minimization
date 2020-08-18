//
//  matrix.h
//  Define Class for Vector & Matrix
//
//  Created by Yoshi Miyazaki on 2014/10/13.
//

#ifndef Gibbs_matrix_h
#define Gibbs_matrix_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include "tensor.h"
using namespace std;

/* return square of a 
template<class T>
inline const T square(const T a){
    return a*a;
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
  }*/

/*----------------------------------------
 1d-vector class
 ---------------------------------------*/
template<class T>
class Vector1d{
private:
    int n;
    T *v;
public:
    Vector1d();
    explicit Vector1d(int nn);
    Vector1d(const T &a, int nn);       // Initialize to constant ca
    Vector1d(const T *a, int nn);       // Initialize to Array
    Vector1d(const Vector1d &copy);      // Copy Constructer
    Vector1d& operator=(const Vector1d &assign);   // Assignment
    Vector1d& operator=(const T &r);               // Assignment w single number
    const bool operator==(const Vector1d<T>& rhs) const;  // bool
    const bool operator!=(const Vector1d<T>& rhs) const{return !(*this == rhs);};
    inline T& operator[](const int i);             // Return i'th element
    inline const T& operator[](const int i) const; // Return i'th element
    inline int  size() const;
    
    /* mathematical operation */
    const T      norm() const;
    const T      maxv() const;
    const T      minv() const;
    const int    maxw() const;
    const int    minw() const;
    const T      sum() const;
    const T      average() const;
    const T      absmaxv() const;
    const T      absminv() const;
    const T      absnon0minv() const;
    const T      absaverage()  const;
    const T      operator*(const Vector1d<T>&) const;    // dot product
    const bool   isnan() const;
    
    const Vector1d<T>  operator+(const Vector1d<T>&);
    const Vector1d<T>  operator+(const Vector1d<T>&) const;
    const Vector1d<T>  operator-(const Vector1d<T>&);
    const Vector1d<T>  operator-(const Vector1d<T>&) const;
    const Vector1d<T>  operator+(const T&);
    const Vector1d<T>  operator+(const T&) const;
    const Vector1d<T>  operator-(const T&);
    const Vector1d<T>  operator-(const T&) const;
    const Vector1d<T>  operator*(const T&);
    const Vector1d<T>  operator*(const T&) const;
    const Vector1d<T>  operator/(const T&);
    const Vector1d<T>  operator/(const T&) const;
    
    Vector1d<T>& operator+=(const Vector1d<T>&);   // add      Vector1d
    Vector1d<T>& operator+=(const T&);             // add         const
    Vector1d<T>& operator-=(const Vector1d<T>&);   // subtract Vector1d
    Vector1d<T>& operator-=(const T&);             // subtract    const
    Vector1d<T>& operator*=(const T&);             // multiply by const
    Vector1d<T>& operator/=(const T&);             // divide   by const

    void resize(int nn);
    void resize(const T& a, int nn);
    void erase (int nn);
    
    tensor1d<T>  to_tensor();
    
    /* destructor */
    ~Vector1d();
};
template<class T>                                            // i'th element
inline T& Vector1d<T>::operator[](const int i){
    return v[i];
}
template<class T>                                            // Const ver.
inline const T& Vector1d<T>::operator[](const int i) const{
    return v[i];
}
template<class T>                                            // Size of Vector
inline int Vector1d<T>::size() const{
    return n;
}

/*----------------------------------------
 2d-vector (matrix) class
 ---------------------------------------*/
template<class T>
class Matrix{ // n row * m column
private:
    int n;
    int m;
    T **v;
public:
    Matrix();
    explicit Matrix(int nn, int mm);
    Matrix(const T &a, int nn, int mm);               // Initialize to constant a
    Matrix(const T *a, int nn, int mm);               // Initialize to Array
    Matrix(const Matrix &copy);                       // Copy Constructer
    Matrix& operator=(const Matrix &copy);            // Assignment
    Matrix& operator=(const T &r);                    // Assignment
    inline T* operator[](const int i);                // Return i'th element
    inline const T* operator[](const int i) const;    // Return i'th element
    inline int nrows() const;
    inline int mcols() const;
    inline T* pointer(){return v[0];}
    
    void resize(int nn, int mm);
    void resize(const T&, int nn, int mm);
    void add_row(Vector1d<double>&);
    void erase_row(int);
    
    Vector1d<T> rowvector(const int j) const;
    Vector1d<T> colvector(const int i) const;
    void        setrowvector(const int, const Vector1d<T>&);
    void        setcolvector(const int, const Vector1d<T>&);

    /* for tensor1d */
    tensor1d<T> rowtensor(const int j) const;
    tensor1d<T> coltensor(const int i) const;
    void        setrowtensor(const int, const tensor1d<T>&);
    void        setcoltensor(const int, const tensor1d<T>&);
    
    /* for arthimetic matrices */
    Matrix<T>   transpose();
    Matrix<T>   lu_decomp();
    void        lu_linear(Vector1d<T>& A);
    Matrix<T>   lu_inverse();
    Matrix<T>&  numeric0(double);
    
    Matrix<T>&  operator+=(const Matrix<T> &B);
    Matrix<T>&  operator-=(const Matrix<T> &B);
    Matrix<T>&  operator*=(const T&);
    Vector1d<T> operator*(Vector1d<T> &A);
    tensor1d<T> operator*(tensor1d<T> &A);
    Matrix<T>   operator*(Matrix<T> &B);

    ~Matrix();
};

template<class T>
inline T* Matrix<T>::operator[](const int i){
    return v[i];
}
template<class T>
inline const T* Matrix<T>::operator[](const int i) const{
    return v[i];
}
template<class T>
inline int Matrix<T>::nrows() const{
    return n;
}
template<class T>
inline int Matrix<T>::mcols() const{
    return m;
}

#include "matrix_s.h"

#endif
