//
//  matrix.cpp
//  Define Class for Vector & Matrix
//
//  Created by Yoshi Miyazaki on 2015/04/11.
//

#include "matrix.h"
/*----------------------------------------
 Vector Types Constructers
 ---------------------------------------*/
template<class T>
Vector1d<T>::Vector1d(){
    n = 0;
    v = 0;
}
template<class T>
Vector1d<T>::Vector1d(int nn){
    n = nn;
    v = new T[n];
}
template<class T>
Vector1d<T>::Vector1d(const T& a, int nn){
    n = nn;
    v = new T[nn];
    for (int i=0; i<nn; i++){
        v[i] = a;
    }
}
template<class T>
Vector1d<T>::Vector1d(const T* a, int nn){
    n = nn;
    v = new T[n];
    for (int i=0; i<nn; i++){
        v[i] = *a++;
    }
}
template<class T>
Vector1d<T>::Vector1d(const Vector1d<T> &copy){
    n = copy.n;
    v = new T[n];
    for (int i=0; i<n; i++){
        v[i] = copy[i];
    }
}
/*----------------------------------------
 Operater
 ---------------------------------------*/
template<class T>
Vector1d<T>& Vector1d<T>::operator=(const Vector1d<T> &copy){
    if (this != &copy){
        if (n != copy.n){
            if (v != 0) delete[] v;
            n = copy.n;
            v = new T[n];
        }
        for (int i=0; i<n; i++){
            v[i] = copy[i];
        }
    }
    return *this;
}
template<class T>
Vector1d<T>& Vector1d<T>::operator=(const T &a){
    for (int i=0; i<n; i++){
        v[i] = a;
    }
    return *this;
}
template<class T>
const bool Vector1d<T>::operator==(const Vector1d<T>& rhs) const{
    if (n != rhs.n){
        return 0;
    }
    else{
        bool b = 1;
        for (int i=0; i<n; i++){
            if (v[i] != rhs[i]){
                b = 0;
                break; 
            }
        }
        return b;
    }
}
template<class T>
void Vector1d<T>::resize(int nn){
    if (n != nn){
        if (v != 0){
            delete[] v;
        }
        n = nn;
        v = new T[n];
    }
}
template<class T>
void Vector1d<T>::resize(const T& a, int nn){
    T *copy = new T[n];
    for (int i=0; i<n; i++){ copy[i] = v[i]; }
    
    int n_old = n;
    if (n != nn){
        if (v != 0){ delete[] v; }
        n = nn;
        v = new T[n];
    }
    for (int i=0; i<n_old; i++){ v[i] = copy[i];}
    for (int i=n_old; i<n; i++){ v[i] = a;      }
    
    if (copy != 0){ delete[] copy; }
}
template<class T>
void Vector1d<T>::erase(int ir){
    if (ir < 0 || n <= ir){ return; } /* if index is outside the range */
    
    T *copy = new T[n];
    for (int i=0; i<n; i++){ copy[i] = v[i]; }
    
    if (v != 0){ delete[] v; }
    n--;  v = new T[n];
        
    for (int i=0;  i<ir; i++){ v[i] = copy[i];   }
    for (int i=ir; i<n;  i++){ v[i] = copy[i+1]; }
    
    if (copy != 0){ delete[] copy; }
}

/*----------------------------------------
 Mathematical Operater
 ---------------------------------------*/
template<class T>
const T Vector1d<T>::norm() const{
    T norm = 0;
    for (int i=0; i<n; i++){
        norm += v[i]*v[i];
    }
    return sqrt(norm);
}
template<class T>
const T Vector1d<T>::maxv() const{
    T maxv = v[0];
    for (int i=1; i<n; i++){
	if (maxv < v[i]){maxv = v[i];}
    }
    return maxv;
}
template<class T>
const T Vector1d<T>::minv() const{
    T minv = v[0];
    for (int i=1; i<n; i++){
	if (minv > v[i]){minv = v[i];}
    }
    return minv;
}
template<class T>
const int Vector1d<T>::maxw() const{
    T maxv = v[0]; int maxw = 0;
    for (int i=1; i<n; i++){
	if (maxv < v[i]){maxv = v[i]; maxw = i;}
    }
    return maxw;
}
template<class T>
const int Vector1d<T>::minw() const{
    T minv = v[0]; int minw = 0;
    for (int i=1; i<n; i++){
	if (minv > v[i]){minv = v[i]; minw = i;}
    }
    return minw;
}
template<class T>
const T Vector1d<T>::sum() const{
    T tot = 0;
    for (int i=0; i<n; i++){ tot += v[i]; }
    return tot;
}
template<class T>
const T Vector1d<T>::average() const{
    T ave = 0;
    for (int i=0; i<n; i++){
	ave += v[i];
    }
    return ave/double(n);
}
template<class T> /* maximum of abs(v[i]) */
const T Vector1d<T>::absmaxv() const{
    T maxv = abs(v[0]);
    for (int i=1; i<n; i++){
	if (maxv < abs(v[i])){maxv = abs(v[i]);}
    }
    return maxv;
}
template<class T> /* minimum of abs(v[i]) */
const T Vector1d<T>::absminv() const{
    T minv = abs(v[0]);
    for (int i=1; i<n; i++){
	if (minv > abs(v[i])){minv = abs(v[i]);}
    }
    return minv;
}
template<class T> /* minimum of abs(v[i]) */
const T Vector1d<T>::absnon0minv() const{
    T minv = absmaxv();
    for (int i=0; i<n; i++){
      if ((minv > abs(v[i])) && (v[i] != 0)){minv = abs(v[i]);}
    }
    return minv;
}
template<class T> /* average of abs(v[i]) */
const T Vector1d<T>::absaverage() const{
    T ave = 0;
    for (int i=0; i<n; i++){
	ave += (v[i]>0 ? v[i] : -1.0*v[i]);
    }
    return ave/double(n);
}
template<class T> /* dot product */
const T Vector1d<T>::operator*(const Vector1d<T>& A) const{
    int nA;    nA = A.size();
    T   dotp = 0;

    if (nA != n){
        cout << "size of vectors don't match (*). Revise your input." << endl;
        exit(7);
    }
    else{
        for (int i=0; i<n; i++){
            dotp += v[i]*A[i];
        }
        return dotp;
    }
}
template<class T>
const bool Vector1d<T>::isnan() const{
    bool isNAN = false;
    for (int i=0; i<n; i++){
	T current = v[i];
	if(std::isnan(current)){ isNAN = true;  break; }
    }
    return isNAN;
}
template<class T>
const Vector1d<T> Vector1d<T>::operator+(const Vector1d<T>& A){
    int nA = A.size();
    if (nA != n){
        cout << "size of vectors don't match (+). Revise your input." << endl;
        exit(7);
    }
    else{
        Vector1d<double> sum(n);
        for (int i=0; i<n; i++){ sum[i] = v[i] + A[i]; }
        return sum;
    }
}
template<class T>
const Vector1d<T> Vector1d<T>::operator+(const Vector1d<T>& A) const{
    int nA = A.size();
    if (nA != n){
        cout << "size of vectors don't match (+). Revise your input." << endl;
        exit(7);
    }
    else{
        Vector1d<double> sum(n);
        for (int i=0; i<n; i++){ sum[i] = v[i] + A[i]; }
        return sum;
    }
}
template<class T>
const Vector1d<T> Vector1d<T>::operator-(const Vector1d<T>& A){
    int nA = A.size();
    if (nA != n){
        cout << "size of vectors don't match (-). Revise your input." << endl;
        exit(7);
    }
    else{
        Vector1d<double> sum(n);
        for (int i=0; i<n; i++){ sum[i] = v[i] - A[i]; }
        return sum;
    }
}
template<class T>
const Vector1d<T> Vector1d<T>::operator-(const Vector1d<T>& A) const{
    int nA = A.size();
    if (nA != n){
        cout << "size of vectors don't match (-). Revise your input." << endl;
        exit(7);
    }
    else{
        Vector1d<double> sum(n);
        for (int i=0; i<n; i++){ sum[i] = v[i] - A[i]; }
        return sum;
    }
}
template<class T>
const Vector1d<T> Vector1d<T>::operator+(const T& A){
    Vector1d<double> sum(n);
    for (int i=0; i<n; i++){
        sum[i] = v[i] + A;
    }
    return sum;
}
template<class T>
const Vector1d<T> Vector1d<T>::operator+(const T& A) const{
    Vector1d<double> sum(n);
    for (int i=0; i<n; i++){
        sum[i] = v[i] + A;
    }
    return sum;
}
template<class T>
const Vector1d<T> Vector1d<T>::operator-(const T& A){
    Vector1d<double> sum(n);
    for (int i=0; i<n; i++){
        sum[i] = v[i] - A;
    }
    return sum;
}
template<class T>
const Vector1d<T> Vector1d<T>::operator-(const T& A) const{
    Vector1d<double> sum(n);
    for (int i=0; i<n; i++){
        sum[i] = v[i] - A;
    }
    return sum;
}
template<class T>
const Vector1d<T> Vector1d<T>::operator*(const T& A){
    Vector1d<double> product(n);
    for (int i=0; i<n; i++){
        product[i] = v[i] * A;
    }
    return product;
}
template<class T>
const Vector1d<T> Vector1d<T>::operator*(const T& A) const{
    Vector1d<double> product(n);
    for (int i=0; i<n; i++){
        product[i] = v[i] * A;
    }
    return product;
}
template<class T>
const Vector1d<T> Vector1d<T>::operator/(const T& A){
    Vector1d<double> quotient(n);
    for (int i=0; i<n; i++){
        quotient[i] = v[i] / A;
    }
    return quotient;
}
template<class T>
const Vector1d<T> Vector1d<T>::operator/(const T& A) const{
    Vector1d<double> quotient(n);
    for (int i=0; i<n; i++){
        quotient[i] = v[i] / A;
    }
    return quotient;
}
template<class T>
Vector1d<T>& Vector1d<T>::operator+=(const Vector1d<T>& A){
    int nA;
    nA = A.size();
    
    if (nA != n){
        cout << "size of vectors don't match (+=). Revise your input." << endl;
        exit(7);
    }
    else{
        for (int i=0; i<n; i++){
            v[i] += A[i];
        }
        return *this;
    }
}
template<class T>
Vector1d<T>& Vector1d<T>::operator+=(const T& a){
    for (int i=0; i<n; i++){
        v[i] += a;
    }
    return *this;
}
template<class T>
Vector1d<T>& Vector1d<T>::operator-=(const Vector1d<T>& A){
    int nA;
    nA = A.size();
    
    if (nA != n){
        cout << "size of vectors don't match (-=). Revise your input." << endl;
        exit(7);
    }
    else{
        for (int i=0; i<n; i++){
            v[i] -= A[i];
        }
        return *this;
    }
}
template<class T>
Vector1d<T>& Vector1d<T>::operator-=(const T& a){
    for (int i=0; i<n; i++){
        v[i] -= a;
    }
    return *this;
}
template<class T>
Vector1d<T>& Vector1d<T>::operator*=(const T& a){
    for (int i=0; i<n; i++){
        v[i] *= a;
    }
    return *this;
}
template<class T>
Vector1d<T>& Vector1d<T>::operator/=(const T& a){
    for (int i=0; i<n; i++){
        v[i] /= a;
    }
    return *this;
}
template<class T>
tensor1d<T> Vector1d<T>::to_tensor(){
    tensor1d<T>  conv(n);
    int i=0;
    for (auto it=conv.begin(); it!=conv.end(); it++){
        *it = v[i];  i++;
    }
    return conv;
}

/*----------------------------------------
 Destructers
 ---------------------------------------*/
template<class T>
Vector1d<T>::~Vector1d<T>(){
    if (v != 0){
        delete[] (v);
    }
}


/*----------------------------------------
 Matrix Types Constructers
 ---------------------------------------*/
template<class T>
Matrix<T>::Matrix(){
    n = 0;    m = 0;
    v = 0;
}
template<class T>
Matrix<T>::Matrix(int nn, int mm){
    n = nn;    m = mm;
    v = new T*[n];
    v[0] = new T[m*n];
    for (int i=1; i<n; i++){
        v[i] = v[i-1] + m;
    }
}
template<class T>
Matrix<T>::Matrix(const T &a, int nn, int mm){
    n = nn;    m = mm;
    v = new T*[n];
    v[0] = new T[m*n];
    for (int i=1; i<n; i++){
        v[i] = v[i-1] + m;
    }
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            v[i][j] = a;
        }
    }
}
template<class T>
Matrix<T>::Matrix(const T *a, int nn, int mm){
    n = nn;    m = mm;
    v = new T*[n];
    v[0] = new T[m*n];
    for (int i=1; i<n; i++){
        v[i] = v[i-1] + m;
    }
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            v[i][j] = *a++;
        }
    }
}
template<class T>
Matrix<T>::Matrix(const Matrix &copy){
    n = copy.n; m = copy.m;
    v = new T*[n];
    v[0] = new T[m*n];
    for (int i=1; i<n; i++){
        v[i] = v[i-1] + m;
    }
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            v[i][j] = copy[i][j];
        }
    }
}

/*----------------------------------------
 Operater
 ---------------------------------------*/
template<class T>
Matrix<T>& Matrix<T>:: operator=(const Matrix<T> &copy){
    if (this != &copy){
        if (n != copy.n || m != copy.m){
            if (v != 0){
                delete v[0];
                delete v;
            }
            n = copy.n;
            m = copy.m;
            v = new T*[n];
            v[0] = new T[n*m];
        }
        for (int i=1; i<n; i++){
            v[i] = v[i-1] + m;
        }
        for (int i=0; i<n; i++){
            for (int j=0; j<m; j++){
                v[i][j] = copy[i][j];
            }
        }
    }
    return *this;
}
template<class T>
Matrix<T>& Matrix<T>:: operator=(const T &r){
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            v[i][j] = r;
        }
    }
    return *this;
}
template<class T>
void Matrix<T>::resize(int nn, int mm){
    if (n != nn || m != mm){
        if (v != 0){
            delete v[0];
            delete v;
        }
        n = nn;
        m = mm;
        v = new T*[n];
        v[0] = new T[n*m];
    }
    for (int i=1; i<n; i++){
        v[i] = v[i-1] + m;
    }
}
template<class T>
void Matrix<T>::resize(const T& a, int nn, int mm){
    if (n != nn || m != mm){
        if (v != 0){
            delete v[0];
            delete v;
        }
        n = nn;
        m = mm;
        v = new T*[n];
        v[0] = new T[n*m];
    }
    for (int i=1; i<n; i++){
        v[i] = v[i-1] + m;
    }
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            v[i][j] = a;
        }
    }
}
template<class T>
void Matrix<T>::add_row(Vector1d<double>& add){
    if (m != add.size()){
        if (m > 0){
            cout << "matrix_s.h: add_row() - vector size unmatch. m = " << m;
            cout << " , add.size() = " << add.size() << endl;
            exit(1);
        } else {
            resize(1,add.size());
            for (int j=0; j<m; j++){ v[0][j] = add[j]; }
            // cout << "row = " << nrows() << " , col = " << mcols() << endl;
            return;
        }
    }
    
    /* copy data to tmp */
    T** tmp = new T*[n];
    tmp[0] = new T[m*n];
    for (int i=1; i<n; i++){ tmp[i] = tmp[i-1] + m; }
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){ tmp[i][j] = v[i][j]; }
    }
    
    /* create new v */
    if (v != 0){
        if (m != 0){ delete[] v[0]; }
        delete[] v;
    }
    n++;
    v = new T*[n];
    v[0] = new T[m*n];
    
    /* copy data */
    for (int i=1; i<n; i++){ v[i] = v[i-1] + m; }
    for (int i=0; i<(n-1); i++){
        for (int j=0; j<m; j++){ v[i][j] = tmp[i][j]; }
    }
    for (int j=0; j<m; j++){ v[n-1][j] = add[j]; }
    
    delete[] tmp[0];
    delete[] tmp;
}
template<class T>
void Matrix<T>::erase_row(int ir){
    if (n == 0){ return; }
    
    /* copy data to tmp */
    T** tmp = new T*[n];
    tmp[0] = new T[m*n];
    for (int i=1; i<n; i++){ tmp[i] = tmp[i-1] + m; }
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){ tmp[i][j] = v[i][j]; }
    }
    
    /* create new v */
    if (v != 0){
        if (m != 0){ delete[] v[0]; }
        delete[] v;
    }
    n--;
    v = new T*[n];  v[0] = new T[m*n];
    for (int i=1; i<n; i++){ v[i] = v[i-1] + m; }
    
    /* copy data */
    for (int i=0; i<ir; i++){
        for (int j=0; j<m; j++){ v[i][j] = tmp[i][j]; }
    }
    for (int i=ir; i<n; i++){
	for (int j=0; j<m; j++){ v[i][j] = tmp[i+1][j]; }
    }
    
    delete[] tmp[0];
    delete[] tmp;
}

/*----------------------------------------
 Return row & column vector
 ---------------------------------------*/
template<class T>
Vector1d<T> Matrix<T>::colvector(const int j) const{
    Vector1d<T> rowv(n);
    for (int i=0; i<n; i++){ rowv[i] = v[i][j];  }
    return rowv;
}
template<class T>
Vector1d<T> Matrix<T>::rowvector(const int i) const{
    Vector1d<T> colv(m);
    for (int j=0; j<m; j++){ colv[j] = v[i][j]; }
    return colv;
}
template<class T>
void Matrix<T>::setrowvector(const int i, const Vector1d<T>& _v){
    for (int j=0; j<m; j++){ v[i][j] = _v[j]; }
}
template<class T>
void Matrix<T>::setcolvector(const int j, const Vector1d<T>& _v){
    for (int i=0; i<n; i++){ v[i][j] = _v[i]; }
}
template<class T>
tensor1d<T> Matrix<T>::coltensor(const int j) const{
    tensor1d<T> rowv(n);
    for (int i=0; i<n; i++){ rowv[i] = v[i][j];  }
    return rowv;
}
template<class T>
tensor1d<T> Matrix<T>::rowtensor(const int i) const{
    tensor1d<T> colv(m);
    for (int j=0; j<m; j++){ colv[j] = v[i][j]; }
    return colv;
}
template<class T>
void Matrix<T>::setrowtensor(const int i, const tensor1d<T>& _v){
    if (m != (int)_v.size()){
        cout << "error in `setrowvector`: wrontg input tensor size. ";
        cout << m << " <-> " << _v.size() << endl;
    }
    for (int j=0; j<m; j++){ v[i][j] = _v[j]; }
}
template<class T>
void Matrix<T>::setcoltensor(const int j, const tensor1d<T>& _v){
    for (int i=0; i<n; i++){ v[i][j] = _v[i]; }
}

/*----------------------------------------
 Mathematical Operater
 ---------------------------------------*/
template<class T>
Matrix<T> Matrix<T>::transpose(){
    Matrix<T> tran(m,n);  int i,j;
    for (i=0; i<n; i++){
	for (j=0; j<m; j++){
	    tran[j][i] = v[i][j];
	}
    }
    return tran;
}
template<class T>
Matrix<T> Matrix<T>::lu_decomp(){
    if (m != n){
        cout << "unable to calculate the inverse" << endl;
        exit(25);
    }
    Matrix<T> lu(m,m);
    /* LU decomposition */
    for (int i=0; i<m; i++){
        /* calculate l_ij */
        for (int j=i; j<m; j++){
            lu[j][i] = v[j][i];
            for (int k=0; k<i; k++){
                lu[j][i] -= lu[k][i]*lu[j][k];
            } 
        }
        /* calculate u_ij */
        for (int j=i+1; j<m; j++){
            lu[i][j] = v[i][j];
            for (int k=0; k<i; k++){
                lu[i][j] -= lu[k][j]*lu[i][k];
            }
            lu[i][j] /= lu[i][i];
        }
    }
    return lu;
}
template<class T>
void Matrix<T>::lu_linear(Vector1d<T>& A){
    /* calculate solution */
    for (int i=0; i<n; i++){
        for (int k=0; k<i; k++){ A[i] -= v[i][k]*A[k]; }
        A[i] /= v[i][i];
    }
    for (int i=n-1; i>=0; i--){
        for (int k=i+1; k<n; k++){
            A[i] -= v[i][k]*A[k];
        }
    }
}
template<class T>
Matrix<T> Matrix<T>::lu_inverse(){
    /* matrix should already been LU decomposed */
    if (m != n){
        cout << "unable to calculate the inverse" << endl;
        exit(25);
    }
    /* prepare identiy matrix */
    Matrix<T> inv(0.0,m,m);
    for (int i=0; i<m; i++){
        inv[i][i] = 1.0;
    }
    /* calculate inverse */
    for (int j=0; j<m; j++){
        for (int i=0; i<n; i++){
            for (int k=0; k<i; k++){ inv[i][j] -= v[i][k]*inv[k][j]; }
            inv[i][j] /= v[i][i];
        }
        for (int i=n-1; i>=0; i--){
            for (int k=i+1; k<n; k++){
                inv[i][j] -= v[i][k]*inv[k][j];
            }
        }
    }
    return inv;
}
template<class T>
Matrix<T>& Matrix<T>::numeric0(double LIM){
    /* find abs max value in matrix */
    T absmaxv = 0.0;
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            if (abs(v[i][j]) > absmaxv) {absmaxv = abs(v[i][j]);}
        }
    }
    /* drop off all numeric error */
    T eps = absmaxv*LIM*16;
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            if (abs(v[i][j]) < eps && v[i][j] != 0){ v[i][j] = 0; }
        }
    }
    return *this;
}
template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& B){
    int nB = B.nrows();
    int mB = B.mcols();

    if ((nB != n) || (mB != m)){
	cout << "size of matrixes don't match (+=). Revise your input." << endl;
	exit(7);
    }
    else {
	for (int i=0; i<n; i++){
	    for (int j=0; j<m; j++){
		v[i][j] += B[i][j];
	    }
	}
	return *this;
    }
}
template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& B){
    int nB = B.nrows();
    int mB = B.mcols();

    if ((nB != n) || (mB != m)){
	cout << "size of matrixes don't match (-=). Revise your input." << endl;
	exit(7);
    }
    else {
	for (int i=0; i<n; i++){
	    for (int j=0; j<m; j++){
		v[i][j] -= B[i][j];
	    }
	}
	return *this;
    }
}
template<class T>
Matrix<T>& Matrix<T>::operator*=(const T& a){
    for (int i=0; i<n; i++){
	for (int j=0; j<m; j++){
	    v[i][j] *= a;
	}
    }
    return *this;
}
template<class T>
Vector1d<T> Matrix<T>::operator*(Vector1d<T> &A){
    int nA;
    nA = A.size();
    // cout << n << m << nB << mB << endl;
    if (nA != m){
        cout << "size of matrix & vector don't match (*). Revise your input. sizes: " << m << " & " << nA << endl;
        exit(7);
    }
    else{
        Vector1d<T> product(n);
        for (int i=0; i<n; i++){
            product[i] = 0;
            for (int k=0; k<m; k++){
                product[i] += v[i][k]*A[k];
            }
        }
        return product;
    }
}
template<class T>
tensor1d<T> Matrix<T>::operator*(tensor1d<T> &A){
    size_t nA = A.size();
    if ((int)nA != m){
        cout << "size of matrix & vector don't match (*). sizes: " << m << " & " << nA << endl;
        exit(7);
    }
    
    else{
        tensor1d<T> product(n);
        for (int i=0; i<n; i++){
            product[i] = 0;
            for (int k=0; k<m; k++){
                product[i] += v[i][k]*A[k];
            }
        }
        return product;
    }
}
template<class T>
Matrix<T> Matrix<T>::operator*(Matrix<T> &B){
    int nB, mB;
    nB = B.nrows(); mB = B.mcols();
    // cout << n << m << nB << mB << endl;
    if (nB != m){
        cout << "size of matricies don't match (*). Revise. " << nB << " x " << m << endl;
        exit(7);
    }
    else{
        Matrix<T> product(n,mB); int i,j,k;
	// int NUM_THREADS=omp_get_num_procs();
	// omp_set_num_threads(NUM_THREADS);
	// #pragma omp parallel for private(j,k)
        for (i=0; i<n; i++){
            for (j=0; j<mB; j++){
                product[i][j] = 0;
                for (k=0; k<m; k++){
                    product[i][j] += v[i][k]*B[k][j];
                }
            }
        }
        return product;
    }
}

/*----------------------------------------
 Destructers
 ---------------------------------------*/
template<class T>
Matrix<T>::~Matrix<T>(){
    if (v!=0){
        if (m!=0){
            delete[] v[0];
        }
        delete[] v;
    }
}
