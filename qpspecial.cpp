
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

extern "C" int dpotrf_(char *UPLO, int* N, double* A, int* LDA, int* INFO);

// extern "C" int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
int main(){
  char uplo = 'U';
  int n = 3;    
  int LDA = n;
  int info;
  int ret;
  
  std::vector<double> A;
  
  A.push_back(2);
  A.push_back(-1);
  A.push_back(0);
  
  A.push_back(-1);
  A.push_back(2);
  A.push_back(-1);
  
  A.push_back(0);
  A.push_back(-1);
  A.push_back(2);
  
  ret = dpotrf_(&uplo, &n, & *A.begin(), &LDA, &info);
  if (0 == ret){
  std::cout << "solution is:";    
  std::cout << "[" << A[0] << ", " << A[1] << ", " << "]" << std::endl;
  std::cout << "Info = " << info << std::endl; 
  }
  return(0);
}
