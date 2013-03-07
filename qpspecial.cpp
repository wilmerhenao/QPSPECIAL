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
/* ----- local macros  -------------------------------------------------- */

//Index must be defined in column major orderA
#define INDEXCM(i, j)     field_start + i * m + j
#define G(i, j)           G[INDEXCM(i, j)]
#define Q(i, j)           Q[INDEXCM(i, j)]
#define QD(i, j)          QD[INDEXCM(i, j)]
#define MAX(A,B)          ((A) > (B)) ? (A) : (B)
#define ABS(A)            ((A) >= 0) ? (A) : -(A)
// LAPACK function declarations
extern "C" void dgemm_(char *, char *, int*, int*,int*, double*, double*, int*, 
		       double*, int*, double*, double*, int*);
extern "C" void sgemm_(char *, char *, int*, int*,int*, float*, float*, int*, 
		       float*, int*, float*, float*, int*);
extern "C" double dlange_(char*, int*, int*, double*, int*, double* );
extern "C" float slange_(char*, int*, int*, float* int*, float*);
extern "C" int dpotrf_(char *UPLO, int* N, double* A, int* LDA, int* INFO);
extern "C" int spotrf_(char *UPLO, int* N, float* A, int* LDA, int* INFO);
extern "C" int dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);
extern "C" int sgetrs_(char*, int*, int*, float*, int*, int*, float*, int*, int*);
// LAPACK function overloading
void mmul_(char *transA, char *transB, int* M, int* N,int* K, double* alpha, double* A,
	   int* LDA, double* B, int* LDB, double* beta, double* C, int* LDC){
  dgemm_(transA, transB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
}

void mmul_(char *transA, char *transB, int* M, int* N,int* K, float* alpha, float* A,
	   int* LDA, float* B, int* LDB, float* beta, float* C, int* LDC){
  sgemm_(transA, transB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
}

double norm_(char* A, int* B, int* C, double* D, int* E, double* F){
  return(dlange_(A, B, C, D, E, F));
}

float norm_(char* A, int* B, int* C, float* D, int* E, float* F){
  return(slange_(A, B, C, D, E, F));
}

int cholesky_(char *UPLO, int* N, double* A, int* LDA, int* INFO){
  return(dpotrf_(UPLO, N, A, LDA, INFO))
}

int cholesky_(char *UPLO, int* N, float* A, int* LDA, int* INFO){
  return(spotrf_(UPLO, N, A, LDA, INFO))
}

int solve_(char* A, int* B, int* C, double* D, int* E, int* F, double* G, int* H, 
	   int* I){
  return(dgetrs_(A, B, C, D, E, F, G, H, I));
}

int solve_(char* A, int* B, int* C, float* D, int* E, int* F, float* G, int* H, 
	   int* I){
  return(sgetrs_(A, B, C, D, E, F, G, H, I));
}

// T can only be float or less than double though.  No higher prec. allowed.
template<typename T>
class qpclass{
private:
  int m, n, field_start, maxit;
  // G points to the matrix G that will be passed to us. e is a vector of ones.
  T *G, *Q, *QD;
  T eta, delta, mu0, tolmu, tolrs, kmu, nQ, krs, ap, ad, y;
  std::vector<T> x, e, z;
  int info;
public:
  qpclass(int m, int n, int, std::vector<T>, T*);
  ~qpclass();
  dotprod(std::vector<T> a, std::vector<T> b);
  copyQD();
  fillGfrompointer(T*);
  HessianfromG();
  optimization();
  postHessian();
  info = iterativePhase();
};

template<typename T>
qpclass<T>::qpclass(int m0, int n0, int maxit0 = 100, std::vector<T> x0,
		    T* G0)
  :m(m0), n(n0), field_start(0), maxit(maxit0), eta(0.9995), delta(3.0), mu0(0.0), 
   tolmu(1e-5), tolrs(1e-5), kmu(0.0), nQ(0,0), krs(0.0), ap(0.0), ad(0.0), y(0.0), 
   info(11){
  assert(m * n);
  if(0 == m * n){
    std::cerr << "The matrix is empty"<< std::endl;
    exit(EXIT_FAILURE);
  }
  std::fill(e.begin(), e.begin() + n, 1.0 / T(n));
  //initialize vector x
  if(0 == x0.size())
    x = e;
  else
    x = x0;
  fillGfromPointer(G0);
  Q = new T[n * n]; // initialize Q (needs to be zero at the beginning)
  for (int i = 0; i < n; i++){
    for(int j = i; j < n; j++){
      Q(i, j) = Q(j, i) = 0.0 //saving a pass
    }
  }
  z = x;
  mu0 = dotprod(x, z);
  mu0 /= T(n);
  kmu = tolmu * mu0;
}

template<typename T>
qpclass<T>::~qpclass(){ 
  delete [] H;
}

template<typename T>
T qpclass<T>::dotprod(std::vector<T> a, std::vector<T> b){
  T res;
  for(std::vector<T>::iterator ait=a.begin(), bit=b.begin(); ait!=a.end();++ait,++bit){
    res += *ait * (*bit);
  }
  return(res);
}

template<typename T>
fillGfrompointer(T *G0){
  G = G0;
}

template<typename T>
qpclass<T>::HessianfromG(){
  char transA = "T";
  char transB = "N";
  T alpha = 1.0, beta = 0.0;
  mmul_(&transA, &transB, &n, &n, &m, &alpha, G, &m, G, &m, &beta, Q, &n);
}

template<typename T>
qpclass<T>::postHessian(){
  // Norm of the Hessian
  char * tnorm = "I";
  nQ = norm_(tnorm, &n, &n, Q, &n, &n) + 2;
  krs = tolrs * nQ;
}

template<typename T>
qpclass<T>::copyQD(){
  for (int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      QD(i, j) = Q(i, j);
    }
  }
}

template<typename T>
qpclass<T>::iterativePhase(){
  T rs, mu = 0.0, r2 = -1.0, M = 0.0, r5, r6, dy;
  //parameters for LAPACK functions
  char* yTrans = "T";
  char* nTrans = "N";
  char* uplo = "U";
  T alpha = 1.0, beta = 0.0;
  int one = 1, info, ret;

  //parameters for interior point method
  std::vector<T> temp;
  std::fill(temp.begin(); temp.end(); 0.0);
  std::vector<T> r1 = temp, r3 = temp, KT = e;
  std::vector<T> zdx = temp;

  for(int k = 0; k < maxit; k++){
    mmult_(nTrans, nTrans, &n, &one, &n, &alpha, Q, &n, x, &n, &beta, temp, &n);
    std::vector<T>::iterator itemp = temp.begin();
    std::vector<T>::iterator ite = e.begin();
    std::vector<T>::iterator itz = z.begin();
    std::vector<T>::iterator itr3 = r3.begin();
    std::vector<T>::iterator itx = x.begin();
    for(std::vector<T>::iterator itr1 = r1.begin; itr1 != r1.end(); 
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
    
    std::vector<T>::iterator itz = z.begin();
    std::vector<T>::iterator itx = x.begin();
    
    int i = 0;
    copyQD();
    for(std::vector<T>::iterator itzdx = zdx.begin(); itzdx != zdx.end(); 
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
    for(std::vector<T>::iterator it = r4.begin(); it != r4.end(); it++, itr3++)
      *it += (*itr3 / *itx); 

    std::vector<T> r7 = r4; // It needs to be started here because r4 will be destroyed
    // next r4 keeps a temporary solution to a system.  So the original r4 ist kaput
    // This is a really bad practice but saves me a lot of memory
    ret = solve_(yTrans, &n, &one, QD, &n, ipiv, & *r4.begin(), &n, &info);
    r5 = dotprod(KT, r4);
    r6 = r2 + r5;
    dy = -r6 / M;
    ite = e.begin();
    for(std::vector<T>::iterator itr7; itr7 != r7.end(); itr7++, ite++)
      *itr7 += *ite * dy;
    ret = solve_(yTrans, &n, &one, QD, &n, ipiv, & *r7.begin(), &n, &info);
    ret = solve_(nTrans, &n, &one, QD, &n, ipiv, & *r7.begin(), &n, &info);
    dx = r7;
    itr3 = r3.begin();
    itz = z.begin();
    itx = x.begin();
    std::vector<T>::operator itdx;
    std::vector<T>::operator itdz = dz.begin();
    for(itdx = dx.begin(); itdx != dx.end(); itdx++, itr3++, itz++, itx++, itdz++)
      *itdz = (*itr3 - *itz * (*itdx))/(*itx);
    /*
      Determine maximal step possible in the approx. tangent direction here primal step
      size
    */
  }
}

template<typename T>
qpclass<T>::optimization(){
  HessianfromG();
  postHessian();
  iterativePhase();
}
