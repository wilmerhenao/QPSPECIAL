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

extern "C" void dpotrf(char *UPLO, int* N, double* A, int* LDA, int* INFO);



int main(){
  char uplo = 'U';
  int n = 3;    
  int LDA = 2;
  int info = 0;

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

  dpotrf(&uplo, &n, & *A.begin(), &LDA, &info);

  std::cout << "solution is:";    
  std::cout << "[" << A[0] << ", " << A[1] << ", " << "]" << std::endl;
  std::cout << "Info = " << info << std::endl; 

  return(0);
}
