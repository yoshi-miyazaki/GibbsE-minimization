/*
  matrix.cpp
  ... define an original class for 1d and 2d-array
  
  created by Yoshi Miyazaki on 2018/07/01.
*/

#include "tensor.h"
using namespace std;

/*--------------------------------------------------------------
 1d-array constructers
 -------------------------------------------------------------*/
template<class T>
tensor1d<T>::tensor1d(size_t n) : vector<T>(n){
}
template<class T>
tensor1d<T>::tensor1d(const T& a, size_t n) : vector<T>(n){
    for (auto it=begin(); it!=end(); it++){ *it = a; }
}
template<class T>
tensor1d<T>::tensor1d(const T* a, size_t n) : vector<T>(n){
    const T* src = a;
    for (auto it=begin(); it!=end(); it++){
        *it = *a;  a++;
    }
}

/*------------------------------------------------------------*/
/* operators */
template<class T>
tensor1d<T>& tensor1d<T>::operator=(const tensor1d& assign){
    if (size() != assign.size()){ resize(assign.size()); }
    auto src = assign.cbegin();
    for (auto it=begin(); it!=end(); it++){
        *it = *src;  src++;
    }
    return *this;
}
template<class T>
tensor1d<T>& tensor1d<T>::operator=(const T& a){
    for (auto it=begin(); it!=end(); it++){
        *it = a;
    }
    return *this;
}
template<class T>
const bool tensor1d<T>::operator==(const tensor1d& rhs) const{
    if (size() != rhs.size()){ return false; }
    auto src = rhs.cbegin();
    for (auto it=cbegin(); it!=cend(); it++){
        if (*it != *src){ return false; }
        src++;
    }
    return true;
}

/*------------------------------------------------------------*/
/* mathematical operation */
template <class T>
const T tensor1d<T>::norm() const{
    T val = 0;
    for (auto it=cbegin(); it!=cend(); it++){ val += (*it)*(*it); }
    return sqrt(val);
}
template <class T>
const T tensor1d<T>::sum() const{
    T tot = 0;
    for (auto it=cbegin(); it!=cend(); it++){ tot += *it; }
    return tot;
}
template <class T>
const T tensor1d<T>::max() const{
    T val = *(cbegin());
    for (auto it=(cbegin()+1); it!=cend(); it++){ val = maxv(val, *it); }
    return val;
}
template <class T>
const T tensor1d<T>::min() const{
    T val = *(cbegin());
    for (auto it=(cbegin()+1); it!=cend(); it++){ val = minv(val, *it); }
    return val;
}
template<class T>
const T tensor1d<T>::absnon0min() const{
    T val = 1e10;
    for (auto it=cbegin(); it!=cend(); it++){
        if (*it == 0) continue;
        val = minv(val, abs(*it));
    }
    return val;
}
template <class T>
const int tensor1d<T>::argmax() const{
    int id = 0, now = 0;  T val = *(cbegin());
    for (auto it=(cbegin()+1); it!=cend(); it++){
        now++;
        if (val < *it){ id = now; val = *it; }
    }
    return id;
}
template <class T>
const int tensor1d<T>::argmin() const{
    int id = 0, now = 0;  T val = *(cbegin());
    for (auto it=(cbegin()+1); it!=cend(); it++){
        now++;
        if (val > *it){ id = now; val = *it; }
    }
    return id;
}
template <class T>
const bool tensor1d<T>::isnan() const{
    for (auto it=cbegin(); it!=cend(); it++){
        if (std::isnan(*it)){ return true; }
    }
    return false;
}
template <class T>
tensor1d<T> tensor1d<T>::operator+(const tensor1d<T>& A) const{
    if (size() != A.size()){ cout << "err in + (size mismatch)" << endl;   exit(2); }
    
    tensor1d<T> res(size());
    
    auto Asrc = A.cbegin();  auto rsrc = res.begin();
    for (auto it=cbegin(); it!=cend(); it++){
        *rsrc = *it + *Asrc;  Asrc++;  rsrc++;
    }
    return res;
}
template <class T>
tensor1d<T> tensor1d<T>::operator-(const tensor1d<T>& A) const{
    if (size() != A.size()){ cout << "err in - (size mismatch)" << endl;   exit(2); }
    
    tensor1d<T> res(size());
    
    auto Asrc = A.cbegin();  auto rsrc = res.begin();
    for (auto it=cbegin(); it!=cend(); it++){
        *rsrc = *it - *Asrc;  Asrc++;  rsrc++;
    }
    return res;
}
template <class T>
T tensor1d<T>::operator*(const tensor1d<T>& A) const{
    if (size() != A.size()){ cout << "err in * (size mismatch)" << size() <<" / "<< A.size() << endl; exit(2); }
    
    T dot = 0;
    auto Asrc = A.cbegin();
    for (auto it=cbegin(); it!=cend(); it++){
        dot += (*it)*(*Asrc);  Asrc++;
    }
    return dot;
}
template <class T>
tensor1d<T> tensor1d<T>::operator/(const tensor1d<T>& A) const{
    if (size() != A.size()){ cout << "err in / (size mismatch)" << endl;   exit(2); }
    
    tensor1d<T> res(size());
    
    auto Asrc = A.cbegin();   auto rsrc = res.begin();
    for (auto it=cbegin(); it!=cend(); it++){
        *rsrc = *it / *Asrc;   Asrc++;  rsrc++;
    }
    return res;
}
template <class T>
tensor1d<T> tensor1d<T>::operator*(const T& var) const{
    tensor1d<T> res(size());
    
    auto rsrc = res.begin();
    for (auto it=cbegin(); it!=cend(); it++){
        *rsrc = (*it) * var;   rsrc++;
    }
    return res;
}
template <class T>
tensor1d<T> tensor1d<T>::operator/(const T& var) const{
    tensor1d<T> res(size());
    
    auto rsrc = res.begin();
    for (auto it=cbegin(); it!=cend(); it++){
        *rsrc = (*it) / var;   rsrc++;
    }
    return res;
}
template <class T>
tensor1d<T>& tensor1d<T>::operator+=(const tensor1d<T>& A){
    if (size() != A.size()){ cout << "err in += (size mismatch)" << endl;  exit(2); }
    
    auto src = A.cbegin();
    for (auto it=begin(); it!=end(); it++){
        *it += *src;  src++;
    }
    return *this;
}
template <class T>
tensor1d<T>& tensor1d<T>::operator-=(const tensor1d<T>& A){
    if (size() != A.size()){ cout << "err in -= (size mismatch)" << endl;  exit(2); }
    
    auto src = A.cbegin();
    for (auto it=begin(); it!=end(); it++){
        *it -= *src;  src++;
    }
    return *this;
}
template <class T>
tensor1d<T>& tensor1d<T>::operator*=(const T& var){
    for (auto it=begin(); it!=end(); it++){
        *it *= var;
    }
    return *this;
}
template <class T>
tensor1d<T>& tensor1d<T>::operator/=(const T& var){
    for (auto it=begin(); it!=end(); it++){
        *it /= var;
    }
    return *this;
}
