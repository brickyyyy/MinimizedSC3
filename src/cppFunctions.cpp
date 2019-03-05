#include <RcppArmadillo.h>

using namespace arma;

//' Compute Euclidean distance matrix by rows
//' 
//' Used in consmx function
//' 
//' @param x A numeric matrix.
// [[Rcpp::export]]
arma::mat ED1(const arma::mat & x) {
	unsigned int outrows = x.n_rows, i = 0, j = 0;
	double d;
	mat out = zeros<mat>(outrows, outrows);

	for (i = 0; i < outrows - 1; i++) {
		arma::rowvec v1 = x.row(i);
		for (j = i + 1; j < outrows; j++) {
			d = sqrt(sum(pow(v1 - x.row(j), 2.0)));
			out(j, i) = d;
			out(i, j) = d;
		}
	}

	return out;
}

//' Compute Euclidean distance matrix by columns
//' 
//' Used in sc3-funcs.R distance matrix calculation
//' and within the consensus clustering.
//' 
//' @param x A numeric matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix ED2(const Rcpp::NumericMatrix & x) {
	unsigned int outcols = x.ncol(), i = 0, j = 0;
	double d;
	Rcpp::NumericMatrix out(outcols, outcols);

	for (j = 0; j < outcols - 1; j++) {
	    Rcpp::NumericVector v1 = x.column(j);
		for (i = j + 1; i < outcols; i++) {
			d = sqrt(sum(pow(v1 - x.column(i), 2.0)));
			out(i, j) = d;
			out(j, i) = d;
		}
	}

	return out;
}

//' Co-association matrix computation
//' 
//' Computes co-association matrix given cluster labels
//'res /= dat.n_cols;
//' 
//' @param dat a matrix containing clustering solutions in columns
// [[Rcpp::export]]
arma::mat consmx(const arma::mat dat) {

	mat res = dat.n_cols * eye<mat>( dat.n_rows, dat.n_rows );

	int i, j, k;
	for (j = 0; j < dat.n_cols; j++) {
		for (i = 0; i < dat.n_rows; i++) {
			for (k = i + 1; k < dat.n_rows; k++) {
				if (dat(i, j) == dat(k, j)) {
				    res(i, k)++;
					res(k, i)++;
				}
			}
		}
	}
	return res;
}

//' Consensus matrix computation
//' 
//' Computes consensus matrix given a co-association matrix
//' 
//' @param matrix a matrix containing co-association matrix
//' @param k number of clusters
// [[Rcpp::export]]
   arma::mat consensus(const arma::mat matrix, int k){
	int n = matrix.n_cols;
	int c = matrix.n_rows;
	
	mat b = mat(n,c*k,fill::zeros);
	int i , j;
	for(i = 0; i < matrix.n_cols; i++){
		for(j = 0; j < matrix.n_rows ; j++){
	     b(i, matrix(j,i)+k*(j-1)) = 1;
		}
    }
	mat res = b.t() / n*c;
  //'add tolerance at convergence=1e-10.
  return (res);
}




//' Graph Laplacian calculation
//' 
//' Calculate graph Laplacian of a symmetrix matrix
//' 
//' @param A symmetric matrix
//' @export
// [[Rcpp::export]]
arma::mat norm_laplacian(arma::mat A) {
    A = exp(-A/A.max());
    arma::rowvec D_row = pow(sum(A), -0.5);
    A.each_row() %= D_row;
    arma::colvec D_col = conv_to< colvec >::from(D_row);
    A.each_col() %= D_col;
    arma::mat res = eye(A.n_cols, A.n_cols) - A;
    return(res);
}

//' Matrix left-multiplied by its transpose
//' 
//' Given matrix A, the procedure returns A'A.
//' 
//' @param x Numeric matrix.
// [[Rcpp::export]]
arma::mat tmult(arma::mat x) {
    return(x.t()*x);
}

