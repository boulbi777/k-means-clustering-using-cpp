#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;

Col<double> Updatedist(const Mat<double> & centers, const Mat<double> x, const int i) {
  Col<double> m_dist = zeros(centers.n_rows);
  for (unsigned int k=0; k<centers.n_rows;k++){
    m_dist(k) = sum(pow(centers.row(k) -  x.row(i), 2));
  }
  return m_dist;
}

double UpdateCriterion(const Mat<double> & centers, const Mat<double> x){
  uword  index;
  double criterion=0;
  Col<double> dist=zeros(centers.n_rows);
  for (unsigned int i=0; i<x.n_rows; i++){
    dist = Updatedist(centers, x, i);
    criterion += dist.min(index);
  }
  return criterion;
}



Mat<double> UpdatePartitions(const Mat<double> & centers, const Mat<double> x){
  Mat<double> zactu = zeros(x.n_rows, centers.n_rows);
  Col<double> dist = zeros( centers.n_rows);
  uword  index;
  double min_val = 0;
  for (unsigned int i=0; i<x.n_rows; i++){
    dist = Updatedist(centers, x, i);
    min_val = dist.min(index);
    zactu(i,index)=1;
  }
  return zactu;
}

Col<double> UpdatePartitionFinal(const Mat<double> & zactu){
  uword  index;
  double min_val = 0;
  Col<double> zactuout=zeros(zactu.n_rows);
  for (unsigned int i=0; i<zactu.n_rows; i++){
    min_val = zactu.row(i).max(index);
    zactuout(i)=index;
  }
  return zactuout + 1;
}

Mat<double> UpdateCenters(const Mat<double> & zactu, const Mat<double> & x){
  Mat<double> centers = trans(zactu) * x ;
  Col<double> effectifs = trans(sum(zactu,0));
  for (unsigned int k=0; k<zactu.n_cols;k++) centers.row(k) = centers.row(k) / effectifs(k);
  return centers;
}

// [[Rcpp::export]]
List singleKmeansC(const arma::mat& x, arma::mat& centers){
  int n = x.n_rows;
  Mat<double> zprec = zeros(n, centers.n_rows), zactu = UpdatePartitions(centers, x);
  while ( sum(sum(zactu % zprec)) <n){
    centers = UpdateCenters(zactu, x);
    zprec = zactu;
    zactu = UpdatePartitions(centers, x);
  }
  double criterion = UpdateCriterion(centers, x);
  return List::create(Named("zactu") = UpdatePartitionFinal(zactu),
                      Named("centers") = centers,
                      Named("criterion") = criterion);
}
