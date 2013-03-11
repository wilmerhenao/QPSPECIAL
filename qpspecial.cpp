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
#include "qpobject.hpp"

int main(){
  
  // This is just here in order to test the stuff from qpobject.hpp and lapackstuff.hpp
  std::cout << "QP special optimizer " << std::endl;
  double* G = new double[9];
  G[1] = 1.0; G[2] = 2.0; G[3] = 3.0; G[0] = 4.0; G[4] = 3.2; G[5] = 21.0; G[6] = 2.21;
  G[7] = 0.0; G[8] = 1.122321;
  qpclass<double> * myopt = new qpclass<double>(3, 3, G, 100);
  myopt->optimization();
  double* solution = new double[3];
  
  myopt->fetchSolution(solution);
  for(int i = 0; i < 3; i++)
    std::cout << solution[i] << std::endl;
  return(0);
}

