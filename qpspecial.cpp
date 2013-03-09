//#define NDEBUG //Uncomment if you are not debugging (FASTER)
/*
  This function sets up the solution of the QP problem
  min             q(x)   = || G*x ||_2^2 = x'*(G'*G)*x
  s.t.            sum(x) = 1
  .                   x >= 0
  The smallest vector under the metric induced by (G' * G) 
  Or the smallest vector in the convex hull of the columns of G under euclidean      
  distance
  Original matlab solution of this problem by Anders Skajaa
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
/* ----- local macros  -------------------------------------------------- */

//Index must be defined in column major orderA
#define INDEXCM(i, j)     field_start + i * m + j
#define G(i, j)           G[INDEXCM(i, j)]
#define Q(i, j)           Q[INDEXCM(i, j)]
#define QD(i, j)          QD[INDEXCM(i, j)]
#define MAX(A,B)          ((A) > (B)) ? (A) : (B)
#define MIN(A,B)          ((A) > (B)) ? (B) : (A)
#define ABS(A)            ((A) >= 0) ? (A) : -(A)
// LAPACK function declarations
extern "C" void dgemm_(char *, char *, int*, int*,int*, double*, double*, int*, 
		       double*, int*, double*, double*, int*);
extern "C" void sgemm_(char *, char *, int*, int*,int*, float*, float*, int*, 
		       float*, int*, float*, float*, int*);
extern "C" double dlange_(char*, int*, int*, double*, int*, double* );
extern "C" float slange_(char*, int*, int*, float*, int*, float*);
extern "C" int dpotrf_(char *UPLO, int* N, double* A, int* LDA, int* INFO);
extern "C" int spotrf_(char *UPLO, int* N, float* A, int* LDA, int* INFO);
extern "C" int dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);
extern "C" int sgetrs_(char*, int*, int*, float*, int*, int*, float*, int*, int*);
// LAPACK function overloading
void mmul_(char transA, char transB, int M, int N,int K, double alpha, double*& A,
	   int LDA, std::vector<double>& B, int LDB, double beta, 
	   std::vector<double>& C, int LDC){
  double* pB = &*B.begin();
  double* pC = &*C.begin();
  dgemm_(&transA, &transB, &M, &N, &K, &alpha, A, &LDA, pB, &LDB, &beta,
	 pC, &LDC);
}

void mmul_(char transA, char transB, int M, int N,int K, float alpha, float* A,
	   int LDA, std::vector<float>& B, int LDB, float beta, 
	   std::vector<float>& C, int LDC){
  sgemm_(&transA, &transB, &M, &N, &K, &alpha, A, &LDA, & *B.begin(), &LDB, &beta, 
	 & *C.begin(), &LDC);
}

void mmul_(char transA, char transB, int M, int N,int K, double alpha, double*& A,
	   int LDA, double* B, int LDB, double beta, 
	   double* C, int LDC){
  dgemm_(&transA, &transB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}

void mmul_(char transA, char transB, int M, int N,int K, float alpha, float*& A,
	   int LDA, float* B, int LDB, float beta, 
	   float* C, int LDC){
  sgemm_(&transA, &transB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}

double norm_(char* A, int* B, int* C, double*& D, int* E, double* F){
  return(dlange_(A, B, C, D, E, F));
}

float norm_(char* A, int* B, int* C, float*& D, int* E, float* F){
  return(slange_(A, B, C, D, E, F));
}

int cholesky_(char &UPLO, int* N, double*& A, int* LDA, int* INFO){
  return(dpotrf_(&UPLO, N, A, LDA, INFO));
}

int cholesky_(char &UPLO, int* N, float*& A, int* LDA, int* INFO){
  return(spotrf_(&UPLO, N, A, LDA, INFO));
}

int solve_(char& A, int* B, int* C, double*& D, int* E, int*& F, double* G, int* H, 
	   int* I){
  return(dgetrs_(&A, B, C, D, E, F, G, H, I));
}

int solve_(char& A, int* B, int* C, float*& D, int* E, int*& F, float* G, int* H, 
	   int* I){
  return(sgetrs_(&A, B, C, D, E, F, G, H, I));
}

/*
template<typename T>
void vecCopy(std::vector<T>& A, std::vector<T>& B){
  // This function will create a deep copy of A in B
  typename std::vector<T>::iterator it;
  for(it = A.begin(); it != A.end(); it++){
    std::cout << "copying " << *it << std::endl;
    B.push_back((T)*it);
  }
  }*/

// T can only be float or less than double though.  No higher prec. allowed.
template<typename T>
class qpclass{
private:
  int m, n, maxit, field_start;
  T eta, delta, mu0, tolmu, tolrs, kmu, nQ, krs, ap, ad, y, q;
  // G points to the matrix G that will be passed to us. e is a vector of ones.
  int info;
  T *G, *Q, *QD;
  std::vector<T> x, e, z, d;
public:
  qpclass(int, int, T*&, int);
  ~qpclass();
  T dotprod(std::vector<T> a, std::vector<T> b);
  void copyQD();
  void fillGfromPointer(T*);
  void HessianfromG();
  void optimization();
  void postHessian();
  int iterativePhase();
  void printSolution();
};
template<typename T>
qpclass<T>::qpclass(int m0, int n0, T*& G0, int maxit0){
  m = m0;
  n = n0;
  maxit = maxit0;
  field_start = 0;
  eta = 0.9995;
  delta = 3.0;
  mu0 = 0.0;
  tolmu = 1e-5;
  tolrs = 1e-5;
  kmu = 0.0;
  nQ = 0.0;
  krs = 0.0;
  ap = 0.0;
  ad = 0.0;
  y = 0.0;
  q = 0.0;
  info = 11;
  assert(m * n);
  if(0 == m * n){
    std::cerr << "The matrix is empty"<< std::endl;
    exit(EXIT_FAILURE);
  }
  for(int i = 0; i < n; i++)
    e.push_back( 1.0/(T)n );
  x = e;

  typename std::vector<T>::iterator ite;
  for(ite = x.begin(); ite != x.end(); ite++){
    std::cout << *ite << " " << std::endl;
    std::cout << &(*ite) << " " << &*e.begin() << std::endl;
  }


  fillGfromPointer(G0);
  Q = new T[n * n]; // initialize Q (needs to be zero at the beginning)
  for (int i = 0; i < n; i++){
    for(int j = i; j < n; j++){
      Q(i, j) = Q(j, i) = 0.0; //saving a pass
    }
  }
  z = x;
  mu0 = dotprod(x, z);
  mu0 /= T(n);
  kmu = tolmu * mu0;
  d = x;
  m = m0;
}

template<typename T>
qpclass<T>::~qpclass(){ 
  delete [] Q;
}

template<typename T>
T qpclass<T>::dotprod(std::vector<T> a, std::vector<T> b){
  T res = 0.0;
  typename std::vector<T>::iterator bit = b.begin();
  for(typename std::vector<T>::iterator ait = a.begin(); ait!=a.end(); ++ait,++bit){
    res += *ait * (*bit);
  }
  return(res);
}

template<typename T>
void qpclass<T>::fillGfromPointer(T *G0){
  G = G0;
}

template<typename T>
void qpclass<T>::HessianfromG(){
  char transA = 'T';
  char transB = 'N';
  T alpha = 1.0, beta = 0.0;
  mmul_(transA, transB, n, n, m, alpha, G, m, G, m, beta, Q, n);
}

template<typename T>
void qpclass<T>::postHessian(){
  // Norm of the Hessian
  char tnorm = 'I';
  T *WORK = new T[n];
  nQ = norm_(&tnorm, &n, &n, Q, &n, WORK) + 2;
  krs = tolrs * nQ;
}

template<typename T>
void qpclass<T>::copyQD(){
  for (int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      QD(i, j) = Q(i, j);
    }
  }
}

template<typename T>
int qpclass<T>::iterativePhase(){
  T rs = 0.0, mu = 0.0, r2 = -1.0, M = 0.0, r5, r6, dy;
  //parameters for LAPACK functions
  char yTrans = 'T';
  char nTrans = 'N';
  char uplo = 'U';
  T alpha = 1.0, beta = 0.0, muaff = 0.0, sig = 0.0;
  int one = 1, ret;
  
  //parameters for interior point method
  std::vector<T> temp;
  std::fill(temp.begin(), temp.end(), 0.0);
  std::vector<T> r1 = temp, r3 = temp, KT = e;
  std::vector<T> zdx = temp, p = temp;
  int k;
  typename std::vector<T>::iterator itx = x.begin();
  for(k = 0; k < maxit; k++){
    mmul_(nTrans, nTrans, n, one, n, alpha, Q, n, x, n, beta, temp, n);
    typename std::vector<T>::iterator itemp = temp.begin();
    typename std::vector<T>::iterator ite = e.begin();
    typename std::vector<T>::iterator itz = z.begin();
    typename std::vector<T>::iterator itr3 = r3.begin();
    typename std::vector<T>::iterator itp = p.begin();
    typename std::vector<T>::iterator itr1;
    typename std::vector<T>::iterator itr4;
    typename std::vector<T>::iterator itzdx;
    typename std::vector<T>::iterator itr7;
    typename std::vector<T>::iterator itdx;
    typename std::vector<T>::iterator itdz;
    for(itr1 = r1.begin(); itr1 != r1.end(); 
	itr1++, itemp++, ite++, itz++, itr3++){
      *itr1 += -(*itemp) + (*ite) * y + (*itz);
      rs = MAX(*itr1, rs);
      r2 += (*itx);
      *itr3 = -((*itx) * (*itz));
      mu -= *itr3;
    }
    rs = MAX(r2, rs); mu /= n;

    if(mu < kmu){
      if(rs < krs){
	return(0); // succesful execution
      }
    }
    
    itz = z.begin();
    itx = x.begin();
    
    int i = 0;
    copyQD();
    for(itzdx = zdx.begin(); itzdx != zdx.end(); 
	itzdx++, itz++, itx++, i++){
      *itzdx = (*itz) / (*itx);
      QD(i, i) += (*itzdx);
    }
    
    // Perform Cholesky factorization on QD
    ret = cholesky_(uplo, &n, QD, &n, &info);
    assert(ret);
    if (0 != info)
      return(2);
    // Solve a system with QD and e (Notice that the solution will be stored in KT)
    int* ipiv = new int[n];
    KT = e; // If I solve directly over e in the next step 'e' will be overwritten!
    ret = solve_(yTrans, &n, &one, QD, &n, ipiv, & *KT.begin(), &n, &info);
    M = dotprod(KT, KT); // might need to make KT members of the class later?

    /* Compute approximate tangent direction using factorization from above */
    std::vector<T> r4 = r1;
    std::vector<T> dx = r1; //r1 just for the size of the vector, not for the values
    std::vector<T> dz = r1;
    itr3 = r3.begin(); itx = x.begin();      
    
    for(itr4 = r4.begin(); itr4 != r4.end(); itr4++, itr3++)
      *itr4 += (*itr3 / *itx); 

    std::vector<T> r7 = r4; // It needs to be started here because r4 will be destroyed
    // next r4 keeps a temporary solution to a system.  So the original r4 ist kaput
    // This is a really bad practice but saves me a lot of memory
    ret = solve_(yTrans, &n, &one, QD, &n, ipiv, & *r4.begin(), &n, &info);
    r5 = dotprod(KT, r4);
    r6 = r2 + r5;
    dy = -r6 / M;
    ite = e.begin();

    for(itr7 = r7.begin(); itr7 != r7.end(); itr7++, ite++)
      *itr7 += *ite * dy;
    ret = solve_(yTrans, &n, &one, QD, &n, ipiv, & *r7.begin(), &n, &info);
    ret = solve_(nTrans, &n, &one, QD, &n, ipiv, & *r7.begin(), &n, &info);
    dx = r7;
    itr3 = r3.begin();
    itz = z.begin();
    itx = x.begin();
    itdz = dz.begin();
    for(itdx = dx.begin(); itdx != dx.end(); itdx++, itr3++, itz++, itx++, itdz++)
      *itdz = (*itr3 - *itz * (*itdx))/(*itx);
    /*
      Determine maximal step possible in the approx. tangent direction here primal step
      size
    */
    itx = x.begin(); itdx = dx.begin(); ap = 1.0;
    itz = z.begin(); itdz = dz.begin(); ad = 1.0;
    for(itp = p.begin(); itp != p.end(); itp++, itx++, itdx++){
      *itp = -(*itx) / *itdx;
      if(*itp > 0 )
	ap = MIN(*itp, ap);
      *itp = -(*itz) / (*itdz);     /* Dual step size */
      if(*itp > 0)  //Using different step sizes in primal and dual improves pfmnce a
	ad = MIN(*itp, ad); //bit
    }
    /* Heuristic for the centering paramater */
    itz = z.begin(); itdx = dx.begin(); itdz = dz.begin();
    for(itx = x.begin(); itx != x.end(); itx++, itz++, itdx++, itdz++){
      muaff += ((*itx) + (ap * (*itdx))) * ((*itz) + (ad * (*itdz)));
    }
    muaff /= (T)n;
    /*
      Compute the new corrected search direction that now includes the appropriate
      amount of centering and mehrotras second order correction term (see r3).  We
      of course reuse the factorization from above
    */
    sig = std::pow((muaff / mu), delta);
    itr3 = r3.begin(); itdx = dx.begin();
    itr4 = r4.begin(); itdz = dz.begin();
    itr1 = r1.begin(); itx = x.begin();
    for(itr3 = r3.begin(); itr3 != r3.end(); itr3++, itdx++, itr4++, itdz++, itr1++,
	  itx++){
      (*itr3) += sig * mu - (*itdx) * (*itdz);
      (*itr4) = (*itr1) + (*itr3) / (*itx);
    }
    r7 = r4;
    ret = solve_(yTrans, &n, &one, QD, &n, ipiv, & *r4.begin(), &n, &info);
    r5 = dotprod(KT, r4);
    r6 = r2 + r5;
    dy = -r6 / M;
    ite = e.begin();
    for( itr7 = r7.begin(); itr7 != r7.end(); itr7++, ite++)
      *itr7 += *ite * dy; //(r7 = r4 + e*dy) -- It is a little cryptic over here!
    ret = solve_(yTrans, &n, &one, QD, &n, ipiv, & *r7.begin(), &n, &info);
    ret = solve_(nTrans, &n, &one, QD, &n, ipiv, & *r7.begin(), &n, &info);
    dx = r7;
    itr3 = r3.begin(); itz = z.begin(); itx = x.begin(); itdz = dz.begin();
    for(itdx = dx.begin(); itdx != dx.end(); itdx++, itr3++, itz++, itx++, itdz++)
      *itdz = (*itr3 - *itz * (*itdx))/(*itx);
    /* Determine maximal step possible in the new direction here primal step size */
    itx = x.begin(); itdx = dx.begin(); ap = 1.0;
    itz = z.begin(); itdz = dz.begin(); ad = 1.0;
    for(itp = p.begin(); itp != p.end(); itp++, itx++, itdx++){
      *itp = -(*itx) / *itdx;
      if(*itp > 0 )
	ap = MIN(*itp, ap);
      *itp = -(*itz) / (*itdz);     /* Dual step size */
      if(*itp > 0)  //Using different step sizes in primal and dual improves pfmnce a
	ad = MIN(*itp, ad); //bit
    }
    /* Update variables primal dual multipliers dual slacks */
    itx = x.begin(); itdx = dx.begin(); itdz = dz.begin();
    for(itz = z.begin(); itz != z.end(); itx++, itz++, itdx++, itdz++){
      (*itx) += eta * ap * (*itdx);    
      y += eta * ad * dy;
      (*itz) += eta * ad * (*itdz);
    }
  }
  if (maxit == k){
    return(1);
  }
  T sumx = 0.0;
  for(itx = x.begin(); itx != x.end(); itx++){
    (*itx) = MAX((*itx), 0.0);
    sumx += (*itx);
  }
  for(itx = x.begin(); itx != x.end(); itx++)
    (*itx) /= sumx;
  d = x;
  mmul_(nTrans, nTrans, m, one, n, alpha, G, n, x, n, beta, d, n);
  q = dotprod(d, d);
  return(0);
}

template<typename T>
void qpclass<T>::printSolution(){
  std::cout << "Solution is " << q << std::endl;
}

template<typename T>
void qpclass<T>::optimization(){
  int ret = 0.0;
  HessianfromG();
  postHessian();
  ret = iterativePhase();
  std::cout << ret << std::endl;
  printSolution();
}

int main(){
  std::cout << "QP special optimizer " << std::endl;
  double* G = new double[4];
  G[1] = 1.0; G[2] = 2.0; G[3] = 3.0; G[0] = 4.0;
  qpclass<double> * myopt = new qpclass<double>(2, 2, G, 100);
  myopt->optimization();
}
