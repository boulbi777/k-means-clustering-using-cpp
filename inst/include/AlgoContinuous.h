#ifndef AlgoContinuous_H
#define AlgoContinuous_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;

class AlgoContinuous{
private:
  unsigned int m_n, m_K;
  uword m_index;
  double m_min_val, m_criterion;
  Mat<double> m_x, m_centers, m_zactu, m_zprec;


public:
  AlgoContinuous(){};
  ~AlgoContinuous(){};
  AlgoContinuous(const arma::mat& x, arma::mat& centers);

  void UpdatePartitions();
  void UpdateCriterion();
  Col<double> UpdatePartitionFinal();
  void UpdateDist();
  void UpdateCenters();
  void Run();
};

#endif

