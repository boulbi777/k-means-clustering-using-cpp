#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List singleKmeansC(const arma::mat& x, arma::mat& centers){
  int K = centers.n_rows, n = x.n_rows;
  uword  index;
  double min_val = 0;
  Col<double> effectifs = zeros(K), dist=zeros(K);
  Mat<double> zprec = zeros(n, K);
  Mat<double> zactu = zeros(n, K);
  for (int i=0; i<n; i++){
    for (int k=0; k<K;k++){
      dist(k) = sum(pow(centers.row(k) - x.row(i), 2));
    }
    min_val = dist.min(index);
    zactu(i,index)=1;
  }
  while ( sum(sum(zactu % zprec)) <n){
    centers = trans(zactu) * x ;
    effectifs = trans(sum(zactu,0));
    for (int k=0; k<K;k++) centers.row(k) = centers.row(k) / effectifs(k);
    zprec = zactu;
    zactu = zactu * 0;
    for (int i=0; i<n; i++){
      for (int k=0; k<K;k++){
        dist(k) = sum(pow(centers.row(k) - x.row(i), 2));
      }
      min_val = dist.min(index);
      zactu(i,index)=1;
    }
  }
  double criterion = 0;
  for (int i=0; i<n; i++){
    criterion += sum(pow(centers.each_row() - x.row(i), 2), 1).min(index);
  }
  Col<double>zactuout = zeros(n);
  for (int i=0; i<n; i++){
    min_val = zactu.row(i).max(index);
    zactuout(i)=index;
  }
  return List::create(Named("zactu") = zactuout + 1,
                      Named("centers") = centers,
                      Named("criterion") = criterion);
}
