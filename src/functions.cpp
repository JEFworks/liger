#include <RcppArmadillo.h>
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
// From Dirk Eddelbuettel (http://gallery.rcpp.org/articles/stl-random-shuffle/)
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// Rcpp implementation of the GSEA resampling procedure
// [[Rcpp::export]]
Rcpp::List gseaRandCore(arma::vec sset, arma::vec eso, int nsamples) {

  std::vector<double> pscores(nsamples);
  std::vector<double> nscores(nsamples);
  int nelem=sset.size();
  
  for(int i=0;i<nsamples;i++) {
    // shuffle the set
    std::random_shuffle(sset.begin(),sset.end(), randWrapper);
    // determine normalizing factors
    double insum=0; double outsum=0;
    for(int j=0;j<nelem;j++) {
      if(sset[j]) {
        insum+=eso[j];
      } else {
        outsum+=eso[j];
      }
    }
    insum/=((double)nelem);
    outsum=(-1.0)*outsum/((double)nelem);
    // calculate score
    double smaxn=0; double smaxp=0; double cs=0;
    for(int j=0;j<nelem;j++) {
      if(sset[j]) {
        cs+=eso[j]/insum;
      } else {
        cs+=eso[j]/outsum;
      }
      if(cs>smaxp) { smaxp=cs; } else if(cs<smaxn) { smaxn=cs; }
    }
    pscores[i]=smaxp; nscores[i]=smaxn;
  }
  return Rcpp::List::create(Rcpp::Named( "p" ) = wrap(pscores),
                            Rcpp::Named( "n" ) = wrap(nscores));
}




// Gsea randomization core method for multiple genesets
// [[Rcpp::export]]
Rcpp::List gseaBulkCore(arma::mat setm, arma::vec eso, int nsamples) {

    int nelem=setm.n_cols;
    int nsets=setm.n_rows;
    
    arma::mat pscores(nsets,nsamples);
    arma::mat nscores(nsets,nsamples);
    
    std::vector<int> colord(nelem);
    for(int i=0;i<nelem;i++) { colord[i]=i; }
    
    arma::vec innv(nsets);
    arma::vec outnv(nsets);
    arma::vec smax(nsets);
    arma::vec smin(nsets);
    arma::vec cs(nsets);
    
    for(int i=0;i<nsamples;i++) {
      // shuffle the order
      std::random_shuffle(colord.begin(),colord.end(), randWrapper);
      
      innv.zeros(); outnv.zeros();
      
      // determine normalizing factors
      for(int j=0;j<nelem;j++) {
        int rj=colord[j];
        for(int k=0;k<nsets;k++) {
          if(setm(k,rj)) {
            innv[k]+=eso[rj];
          } else {
            outnv[k]+=eso[rj];
          }
          
        }
      }
      innv*=(1.0)/((double)nelem);
      outnv*=(-1.0)/((double)nelem);
      // calculate score
      
      smin.zeros(); smax.zeros(); cs.zeros();
      for(int j=0;j<nelem;j++) {
        int rj=colord[j];
        // update cumulative
        for(int k=0;k<nsets;k++) {
          if(setm(k,rj)) {
            cs[k]+=eso[rj]/innv[k];
          } else {
            cs[k]+=eso[rj]/outnv[k];
          }
          if(cs[k]>smax[k]) { smax[k]=cs[k]; } else if(cs[k]<smin[k]) { smin[k]=cs[k]; }
        }
      }
      pscores.col(i)=smax;
      nscores.col(i)=smin;
    }
    return Rcpp::List::create(Rcpp::Named( "p" ) = wrap(pscores),
                              Rcpp::Named( "n" ) = wrap(nscores));
}
  
