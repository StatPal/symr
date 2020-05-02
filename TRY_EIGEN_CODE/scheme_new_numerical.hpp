/*
Precompile header using

g++ scheme_new.hpp -I /usr/include/eigen3 -O3

- Don't do it now. taking huge gch file. 
*/



//#include <RcppEigen.h>

#ifndef MAIN_HEADER
#define MAIN_HEADER

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <random>
#include "meta.h"


#include <cmath> 	//For bessel fn if cpp17
#include <chrono>



// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

//using namespace Rcpp;
using namespace Eigen;

typedef std::function<float(const Vector_eig &x)> my_function_type;
typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::Matrix3f Matrix3f_eig;

#define SQ(x) ((x) * (x))

const int IF_DEBUG = 1;
#define ABS_ERR 1e-10



#define DEBUG_LEVEL_0				//Almost always debug.

#ifdef DEBUG_LEVEL_0
#define Debug0(x) {std::cout << "DEBUG 0: "<< x << "\n";}
#else
#define Debug0(x)
#endif


#define DEBUG_LEVEL_1				//Important ones.

#ifdef DEBUG_LEVEL_1
#define Debug1(x) {std::cout << "DEBUG 1: "<< x << "\n";}
#else
#define Debug1(x)
#endif


//#define DEBUG_LEVEL_2				//All Debugs

#ifdef DEBUG_LEVEL_2
#define Debug2(x) {std::cout << "DEBUG 2: "<< x << "\n";}
#else
#define Debug2(x)
#endif

#define DEBUG_ANOTHER_LEVEL 0


#define DEBUG_LEVEL_times

#ifdef DEBUG_LEVEL_times
#define Debug_time(x) {std::cout << "DEBUG 1: "<< x << "\n";}
#else
#define Debug_time(x)
#endif



/********************************************************
******************** About Eigen ************************
********************************************************/
/*
* .resize may be useful or reshape
* https://gist.github.com/gocarlos/c91237b02c120c6319612e42fa196d77
* x.size() gives ncol*nrow
*/ 



/********************************************************
******************** Bessel functions *******************
********************************************************/

/**
* Modified Bessel Function of First Kind
*/
static double besseli(double x, int alpha){
	double y;
	long double fct1 = 1., sum = 0., tmp = 1.;
	if (x < 15.) {
		y = pow(x/2,2);
		tmp = pow(x/2, alpha)/tgamma(alpha+1);
		sum = tmp;
		for(int k = 1; k < 25; ++k) {
			tmp *= y / (k*(k+alpha));
			sum += tmp; 	//std::cout << "sum = " << sum << "\t\t\tMonitor(" << tmp << ")\n";
			if(fabs(tmp)<ABS_ERR)
				break;
		}
	} else {
		y=(8*x);
		sum = 1., tmp =1.;							//double tmp2 = std::pow(2*alpha, 2);
		int max_k = 10;
		if(alpha>=2)
			max_k = 2*(alpha+(int) sqrt(40*alpha));
		
		for(int k = 1; k < max_k; fct1*=++k) {		// Weird behaviour -- k<50 for nu=0.
			tmp *= ((2*k-1-2*alpha)*(2*k-1+2*alpha))/(k*y);
			sum += tmp; 	//std::cout << "sum = " << sum << "\t\t\tMonitor(" << tmp << ")\n";
			if(fabs(tmp)<ABS_ERR)
				break;
		}
		sum*=exp(x)/sqrt(2*M_PI*x);
	}
	return sum;
}


double our_bessel_I(double x, double nu){
	if(x>=0.0){
	/*
#ifdef (__cplusplus == 201703L)
	return(std::cyl_bessel_i(nu, x));
#else
	return(besseli(x, nu));
#endif*/
		
		return(besseli(x, nu));					//return(R::bessel_i(x, nu, 1));
	} else {
		std::cout << "neg arg in bessel function" << std::endl;
		return -1.0e-50;
	}
}								 //besselI(0.05, 0.0); our_bessel_I(0.05, 0.0)		// value out of range in 'bessel_i'	 -- check


/**
* Ratio of BesselI(x, 1) and BesselI(x, 0)
*/
double ratio_bessel_10(double x){
	if(x<150){
		return (our_bessel_I(x, 1)/our_bessel_I(x, 0));
		
	} else {
		return (1.0 - 0.5/x);
	}
}

/**
* Ratio of BesselI(x, 2) and BesselI(x, 0)
*/
double ratio_bessel_20(double x){
	if(x<150){
		//return (our_bessel_I(x, 2)/our_bessel_I(x, 0));
		return (std::cyl_bessel_i(2, x)/std::cyl_bessel_i(0, x));
	} else {
		return (1.0 - 2/x);		// Take one more term
	}
}

double besselI1_I0(double x)
{ 
  if (x > 3.75) {
    /*
     * Based off of 9.8.2 on p. 378 of Abramowitz and Stegun,
     * currently available at http://www.math.sfu.ca/~cbm/aands/page_378.htm
     * If less precision is required, we can simply throw out
     * some of the v terms, starting with v8 and working down.
     */
    const double tRecip = 3.75 / x;

    const double v8 = 0.00392377 * tRecip;
    const double v7 = (v8 - 0.01647633) * tRecip;
    const double v6 = (0.02635537 + v7) * tRecip;
    const double v5 = (-0.02057706 + v6) * tRecip;
    const double v4 = (0.00916281 + v5) * tRecip;
    const double v3 = (-0.00157565 + v4) * tRecip;
    const double v2 = (0.00225319 + v3) * tRecip;
    const double v1 = 0.39894228 + (0.01328592 + v2) * tRecip;

    /**
     * Based off of 9.8.4 on p. 378 of Abramowitz and Stegun
     */
    const double w8 = -0.00420059 * tRecip;
    const double w7 = (w8 + 0.01787654) * tRecip;
    const double w6 = (w7 - 0.02895312) *tRecip;
    const double w5 = (w6 + 0.02282967) * tRecip;
    const double w4 = (w5 - 0.01031555) * tRecip;
    const double w3 = (w4 + 0.00163801) * tRecip;
    const double w2 = (w3 - 0.00362018) * tRecip;
    const double w1 = 0.39894228 + (w2 - 0.03988024) * tRecip;

    return w1 / v1;
  }
  else if (x == 0.0) {
    return 0.0;
  }
  else {
    /*
     * Based off of 9.8.1 on p. 378 of Abramowitz and Stegun.
     */
    const double t = x/3.75;
    const double tSq = t * t;
    const double v12 = 0.0045813 * tSq;
    const double v10 = (0.0360768 + v12) * tSq;
    const double v8 = (0.2659732 + v10) * tSq;
    const double v6 = (1.2067492 + v8) * tSq;
    const double v4 = (3.0899424 + v6) * tSq;
    const double v2 = 1.0 + (3.5156229 + v4) * tSq;

    /**
     * Based off of 9.8.3 on p. 378 of Abramowitz and Stegun.
     */
    const double w12 = 0.00032411 * tSq;
    const double w10 = (w12 + 0.00301532) * tSq;
    const double w8 = (w10 + 0.02658733) * tSq;
    const double w6 = (w8 + 0.15084934) * tSq;
    const double w4 = (w6 + 0.51498869) * tSq;
    const double w2 = 0.5 + (w4 + 0.87890594) * tSq;
    
    return x * (w2 / v2);
  }

}

double logBesselI0(double x) {
	if (x >= 3.75) {
		/*
		 * Based off of 9.8.2 on p. 378 of Abramowitz and Stegun,
		 * currently available at http://www.math.sfu.ca/~cbm/aands/page_378.htm
		 * If less precision is required, we can simply throw out
		 * some of the v terms, starting with v8 and working down.
		 */
		const double tRecip = 3.75 / x;

		const double v8 = 0.00392377 * tRecip;
		const double v7 = (v8 - 0.01647633) * tRecip;
		const double v6 = (0.02635537 + v7) * tRecip;
		const double v5 = (-0.02057706 + v6) * tRecip;
		const double v4 = (0.00916281 + v5) * tRecip;
		const double v3 = (-0.00157565 + v4) * tRecip;
		const double v2 = (0.00225319 + v3) * tRecip;
		const double v1 = (0.01328592 + v2) * tRecip;
		/*
		 * Is sqrt cheaper than log?
		 * -0.5 log(x) + log(c) = log(c/sqrt(x)) -- inconclusive so far
		 */ 
		return x - 0.5 * log(x) + log(v1 + 0.39894228);
	}
	else if (x == 0.0) {
		return 0.0;
	}
	else {
		/*
		 * Based off of 9.8.1 on p. 378 of Abramowitz and Stegun.
		 */
		 
		const double t = x/3.75;
		const double tSq = t * t;
		const double v12 = 0.0045813 * tSq;
		const double v10 = (0.0360768 + v12) * tSq;
		const double v8 = (0.2659732 + v10) * tSq;
		const double v6 = (1.2067492 + v8) * tSq;
		const double v4 = (3.0899424 + v6) * tSq;
		const double v2 = (3.5156229 + v4) * tSq;
		return log(1.0 + v2);
	}
}



/********************************************************
******************** Other functions ********************
********************************************************/


// [[Rcpp::export]]
Matrix_eig Cov_1(Matrix_eig x) {
	int nRows = x.rows();
	Vector_eig colMean = x.colwise().mean();	// https://eigen.tuxfamily.org/dox/group__TutorialReductionsVisitorsBroadcasting.html
	x.rowwise() -= colMean.transpose();			// or x_cen instead of x
	return x.transpose()*x/(nRows-1);			// or x_cen
}

void show_head(Eigen::MatrixXf W, int n = 10){
	std::cout << "head of the matrix:\n" ;
	for(int i = 0; i < n; ++i){
		std::cout << W.row(i) << "\n";
	}
}

double mean_rice(double nu, double sigma){
	double x = - SQ(nu)/(2*SQ(sigma));
	return sigma * std::sqrt(M_PI/2) * std::exp(x/2)*( (1-x)*our_bessel_I(-x/2, 0) - x * our_bessel_I(-x/2, 1)) ;
}


//[[Rcpp::export]]
Matrix_eig to_matrix(Vector_eig v1, int nrow_in, int ncol_in){
	return(Map<MatrixXf> (v1.data(), nrow_in, ncol_in));
	//https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
	//https://stackoverflow.com/questions/52261389/how-to-convert-an-stdvector-to-a-matrix-in-eigen
	//https://stackoverflow.com/questions/32452739/vector-to-matrix/32475129
}

Vector_eig to_vector(Matrix_eig v1, int is_transpose=0){
	if(!is_transpose){
		return(Map<VectorXf> (v1.data(), v1.rows()*v1.cols()));
	} else {
		return(Map<VectorXf> (v1.transpose().data(), v1.rows()*v1.cols()));
	}
}

// [[Rcpp::export]]
double sp_det_1(SpMat A){
	return MatrixXf(A).determinant();
	//https://stackoverflow.com/questions/15484622/how-to-convert-sparse-matrix-to-dense-matrix-in-eigen/15592295
	//https://stackoverflow.com/questions/13033694/eigen-convert-dense-matrix-to-sparse-one
}

Vector_eig log_vec(Vector_eig x){
	for(int i = 0; i < x.size(); ++i){
		x(i) = std::log(x(i));
	}
	return x;
}

double abs_sum(Vector_eig x){
	double temp = 0.;
	for(int i = 0; i < x.size(); ++i){
		temp += std::fabs(x(i));
	}
	return temp;
}


// [[Rcpp::export]]
double sp_log_det_2(SpMat B){			// Log determinant
	MatrixXf A = MatrixXf(B);
	SelfAdjointEigenSolver<MatrixXf> es(A);				// Marked was in eigen 2
	return log_vec(es.eigenvalues()).sum();
}		// Very small Negative eigenvalue due to numerical error possibly 

double sp_log_det_7(SpMat A){			// Log determinant - LU
	Eigen::SparseLU<Eigen::SparseMatrix<float>, COLAMDOrdering<int>> solver;
	Debug2("Solver initiated!");
	A.makeCompressed();
	Debug2("A compressed!");
	solver.compute(A);
	Debug2("Solver computed!");
	// broken in pieces:
		//solver.analyzePattern(A);
		//solver.factorize(A);
	double temp = solver.logAbsDeterminant();
	Debug2("logAbsDeterminant calculated!");
	return temp;
}



// [[Rcpp::export]]
double log_det_2(Matrix_eig B){			// Log determinant
	SelfAdjointEigenSolver<MatrixXf> es(B);
	return log_vec(es.eigenvalues()).sum();
}

double log_det_3(Matrix3f_eig B){			// Log determinant for fixed size matrix
	SelfAdjointEigenSolver<Matrix3f> es(B);
	return log_vec(es.eigenvalues()).sum();
}

double sp_trace(SpMat A){
	double sum =0;
	for (int k = 0; k < A.outerSize(); ++k) {
		sum += A.coeff(k,k);
	}
	return sum;
}


/********************************************************
******** New Functiones needed for Likelihood ***********
********************************************************/



// [[Rcpp::export]]
SpMat I_n(int n_x){
	SpMat temp(n_x, n_x);				//temp.reserve(n_x);
	temp.setIdentity();
	return(temp);
}

// [[Rcpp::export]]
SpMat J_n(int n_x){				// has determinant 0???
	SpMat temp(n_x, n_x);
	temp.reserve(3*n_x-4);
	for(int i = 0; i < n_x-1; i++){
		temp.insert(i, i) = 2.0;
		temp.insert(i,i+1) = -1;
		temp.insert(i+1,i) = -1;
	}
	temp.coeffRef(0,0) = 1;
	temp.insert(n_x-1,n_x-1) = 1;
	return(temp);
}

Matrix_eig Kron_eig(Matrix_eig m1, Matrix_eig m2){						//https://forum.kde.org/viewtopic.php?f=74&t=50952
	//std::cout << "Kroneker prod started\n";
	Matrix_eig m3(m1.rows()*m2.rows(), m1.cols()*m2.cols());
	
	for (int i = 0; i < m1.cols(); i++) {
		for (int j = 0; j < m1.rows(); j++) {
			m3.block(i*m2.rows(), j*m2.cols(), m2.rows(), m2.cols()) =  m1(i,j)*m2;
		}
	}
	//std::cout << "Kroneker prod done\n";
	return m3;
}

Vector_eig Kron_vec_eig(Vector_eig m1, Vector_eig m2){
	//std::cout << "Kroneker prod started\n";
	Vector_eig m3(m1.size()*m2.size());
	for (int i = 0; i < m1.size(); i++) {
		m3.segment(i*m2.size(), m2.size()) =  m1(i) * m2;
	}
	//std::cout << "Kroneker prod done\n";
	return m3;
}

SpMat Kron_Sparse_eig(SpMat m1, SpMat m2){
	m1.makeCompressed(); m2.makeCompressed();
	int m = m1.rows(), n = m1.cols(), p = m2.rows(), q = m2.cols(), nz1 = m1.nonZeros(), nz2 = m2.nonZeros();
	
	SpMat m3(m*p, n*q);
	m3.reserve(nz1*nz2);
	
	for (int k1=0; k1<m1.outerSize(); ++k1){
		for (int k2=0; k2<m2.outerSize(); ++k2){
			for (SpMat::InnerIterator it1(m1,k1); it1; ++it1){
				for (SpMat::InnerIterator it2(m2,k2); it2; ++it2){
					m3.insert(p * it1.row() + it2.row(), q * it1.col() + it2.col()) = it1.value()*it2.value();
				}
			}
		}
	}
	m3.makeCompressed();
	return m3;
}



// https://stackoverflow.com/questions/30188482/sparse-eigenvalues-using-eigen3-sparse
Vector_eig eigenvals_J_n(int n){
	Matrix_eig A = Matrix_eig(J_n(n));
	Eigen::SelfAdjointEigenSolver<Matrix_eig> es(A);
	return (es.eigenvalues());
}


double sp_log_det_specific(Vector_eig beta, int n_x, int n_y, int n_z, double thres = 0.000001){
	Debug2("Log det calculation started");
	Vector_eig eigens(n_x*n_y*n_z); 
	eigens = beta(0) * Kron_vec_eig(eigenvals_J_n(n_x), Vector_eig::Ones(n_y*n_z)) +
						beta(1)*Kron_vec_eig(Vector_eig::Ones(n_x), Kron_vec_eig(eigenvals_J_n(n_y), Vector_eig::Ones(n_z)) )+
						beta(2) * Kron_vec_eig(Vector_eig::Ones(n_x*n_y), eigenvals_J_n(n_z));
	
	double temp = 0.0;
	for(int i = 0; i < eigens.size(); ++i){
		if(eigens(i)>thres){
			temp += std::log(eigens(i));
		}
	}
	return temp;
}

double sp_log_inv_specific(Vector_eig beta, int n_x, int n_y, int n_z, int k, double thres = 0.000001){
	int n = n_x*n_y*n_z;
	Vector_eig eigens(n); 
	
	Vector_eig D_x = beta(0) * eigenvals_J_n(n_x);
	Vector_eig D_y = beta(1) *eigenvals_J_n(n_y);
	Vector_eig D_z = beta(2) *eigenvals_J_n(n_z);
	
	eigens =  Kron_vec_eig(D_x, Vector_eig::Ones(n_y*n_z)) + Kron_vec_eig(Vector_eig::Ones(n_x), Kron_vec_eig(D_y, Vector_eig::Ones(n_z)) ) +
						Kron_vec_eig(Vector_eig::Ones(n_x*n_y), D_z);
	
	Vector_eig final(n-1);
	if(k == 0){
		final = Kron_vec_eig(D_x, Vector_eig::Ones(n_y*n_z)).segment(0, n-1) ;
	} else if(k == 1){
		final = Kron_vec_eig(Vector_eig::Ones(n_x), Kron_vec_eig(D_y,Vector_eig::Ones(n_z)) ).segment(0, n-1);
	} else if(k == 2){
		final = Kron_vec_eig(Vector_eig::Ones(n_x*n_y), D_z).segment(0, n-1);
	}
	for(int i = 0; i < n-1; ++i){
		final(i) /= eigens(i);
	}
	double temp = (double)final.sum();
	return (temp);
}


// see also https://stackoverflow.com/questions/38839406/eigen-efficient-kronecker-product
// or kroneckerProduct
// [[Rcpp::export]]
SpMat Lambda(Vector_eig beta, int n_x, int n_y, int n_z){
	return(beta(0)*Kron_Sparse_eig(J_n(n_x), I_n(n_y*n_z)) + 
			 beta(1)*Kron_Sparse_eig( Kron_Sparse_eig(I_n(n_x), J_n(n_y)), I_n(n_z)) + 
			 beta(2)*Kron_Sparse_eig(I_n(n_x*n_y), J_n(n_z)));
}

//[[Rcpp::export]]
Eigen::VectorXf Bloch_vec(Eigen::VectorXf W, Eigen::VectorXf TE, Eigen::VectorXf TR){
	int m = TE.size();
	Eigen::VectorXf tmp = Eigen::VectorXf::Zero(m);
	if(W(1)<=0){										// To avoid numerical problem at boundary. // changed to <= 
		for(int j = 0; j < m; ++j) {
			tmp(j) = W(0)*exp(TE(j)*log(W(2)));
		}
	} else {
		for(int j = 0; j < m; ++j) {
			tmp(j) = W(0)*std::exp(TE(j)*std::log(W(2)))*(1-std::exp(TR(j)*std::log(W(1))));
		}
	}
	return tmp;
}

void check_bounds(Matrix_eig &W, Vector_eig lb, Vector_eig ub){
	int n = W.rows();
	int count = 0;
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < 3; ++j){
			Debug2(W(i,j));
			/*try{
				if(W(i, j) < lb(j*n + i)){
					W(i, j) = lb(j*n + i);
					Debug2("Bound Checks: i:" << i << " j:" << j << " W: " << W(i,j));
					count++;
				} else if(W(i, j) > ub(j*n + i)){
					W(i, j) = ub(j*n + i);
					Debug2("Bound Checks: i:" << i << " j:" << j << " W: " << W(i,j));
					count++;
				}
			} catch (...){
				Debug0("i: "<< i << ", j:" << j << "; j*n + i: " << j*n + i);
			}
			*/
		}
	}
	if(count>0){
		Debug1(" "<< count << "boundary cases");
	}
}


int check_nan(Matrix_eig A){
	for(int i = 0; i < A.rows(); ++i){
		for(int j = 0; j < A.cols(); ++j){
			if(std::isnan(A(i, j))){
				Debug0("NAN in location: ("<< i << "," << j<< ")!");
				return 1;
			}
		}
	}
	return 0;
}


int check_nan_vec(Vector_eig A){
	int count = 0;
	for(int i = 0; i < A.size(); ++i){
		if(std::isnan(A(i))){
			Debug0("NAN in location: ("<< i << ")!");
			return 1;
			count++;
		}
	}
	if(count > 0){
		return 1;
	} else{ 
		return 0;
	}
}


// [[Rcpp::export]]
Matrix_eig v_mat(Matrix_eig W, Vector_eig TE, Vector_eig TR){
	int nCol = TE.size();	//m
	int nRow = W.rows();	//n
	Matrix_eig tmp = Matrix_eig::Zero(nRow, nCol);		//check the order
	for(int i = 0; i < nRow; ++i) {
		for(int j = 0; j < nCol; ++j) {
			if(W(i, 1)<=0){								//changed
				tmp(i,j) = W(i,0)*std::exp(TE(j)*std::log(W(i,2)));
			} else {
				tmp(i,j) = W(i,0)*std::exp(TE(j)*std::log(W(i,2)))*(1-std::exp(TR(j)*std::log(W(i,1))));
			}
			if(std::isnan(tmp(i, j))){
				Debug0("NAN in v_mat i:" << i << " j:" << j << " W: " << W(i,0) << ", " << W(i,1) << ", " << W(i,2)  << 
				"\t others:" << std::exp(TE(j)*std::log(W(i,2))) << ", " << (1-std::exp(TR(j)*std::log(W(i,1)))) );
				exit(EXIT_FAILURE);
			}
		}
	}
	return tmp;
}


//[[Rcpp::export]]
Matrix_eig to_W(Vector_eig rho, Vector_eig T_1, Vector_eig T_2){
	Matrix_eig W = Matrix_eig::Zero(rho.size(), 3);
	W.col(0) = rho;				//as<arma::vec>(wrap(rho));
	for(int i = 0; i < rho.size(); ++i){
		W(i, 1) = exp(-1/T_1(i));
		W(i, 2) = exp(-1/T_2(i));
	}
	return(W);
}


Vector_eig to_param_vec(Matrix_eig W, Matrix3f_eig Psi_inv, double beta_x, double beta_y){
	int n = W.rows();		// check
	Vector_eig temp = Vector_eig::Zero(3*n+6+2);

	//int ind =0;
	//for(int j = 0; j < 3; ++j){
	//	for(int i = 0; i < n; ++i){
	//		temp(ind) = W(i, j);
	//		ind++;
	//	}
	//}

	temp.segment(0,3*n)= to_vector(W); // W.resize(n*3, 1) //vectorise( W, 0 );	//Eigen //http://eigen.tuxfamily.org/dox/group__QuickRefPage.html
	temp(3*n+0) = Psi_inv(0,0); temp(3*n+1) = Psi_inv(0,1); temp(3*n+2) = Psi_inv(0,2);
								temp(3*n+3) = Psi_inv(1,1); temp(3*n+4) = Psi_inv(1,2);
															temp(3*n+5) = Psi_inv(2,2);
	temp(3*n+6) = beta_x;
	temp(3*n+7) = beta_y;

	return temp;
}


Vector_eig to_param_vec_grad(Matrix_eig W, Matrix_eig Psi_inv, double beta_x, double beta_y){
	int n = W.rows();		// check
	Vector_eig temp = Vector_eig::Zero(3*n+6+2);
	temp.segment(0,3*n)= to_vector(W); // W.resize(n*3, 1) //vectorise( W, 0 );	//Eigen //http://eigen.tuxfamily.org/dox/group__QuickRefPage.html
	temp.segment(3*n, 6) = Psi_inv;
	temp(3*n+6) = beta_x;
	temp(3*n+7) = beta_y;

	return temp;
}



Matrix3f_eig to_Psi_inv(Vector_eig temp_psi){
	Matrix3f_eig Psi_inv = Matrix3f_eig::Zero(3,3);
	Psi_inv(0,0) = temp_psi(0); Psi_inv(0,1) = temp_psi(1); Psi_inv(0,2) = temp_psi(2);
	Psi_inv(1,0) = temp_psi(1); Psi_inv(1,1) = temp_psi(3); Psi_inv(1,2) = temp_psi(4);
	Psi_inv(2,0) = temp_psi(2); Psi_inv(2,1) = temp_psi(4); Psi_inv(2,2) = temp_psi(5);
	return Psi_inv;
}


Matrix3f_eig to_L_mat(Vector_eig temp_L){
	Matrix3f_eig L = Matrix3f_eig::Zero(3,3);
	L(0,0) = temp_L(0); L(0,1) = temp_L(1); L(0,2) = temp_L(2);
	L(1,0) = temp_L(1); L(1,1) = temp_L(3); L(1,2) = temp_L(4);
	L(2,0) = temp_L(2); L(2,1) = temp_L(4); L(2,2) = temp_L(5);
	return L;
}


Matrix3f_eig from_Cholesky(Matrix3f_eig L){
	return (L*L.transpose());
}

Matrix3f_eig to_Cholesky(Matrix3f_eig A){
	
	Matrix3f_eig L;
	
	L(0,0) = std::sqrt(A(0,0));		L(0,1) = 0; 								L(0,2) = 0;
	L(1,0) = A(1,0)/L(0,0);			L(1,1) = std::sqrt(A(1,1)-SQ(L(1,0)));		L(1,2) = 0;
	L(2,0) = A(2,0)/L(0,0);			L(2,1) = (A(2,1)-L(2,0)*L(1,0))/L(1,1);		L(2,2) = std::sqrt(A(2,2)-SQ(L(2,0))-SQ(L(2,1)));
	
	return L;
}



Matrix_eig to_grad_Cholesky(Vector_eig L){
	
	Matrix_eig D = Matrix_eig::Zero(6, 6);
	
	D(0,0) = 2*L(0);	D(0,1) = L(1);		D(0,2) = L(2);
	D(1,1) = L(0);		D(1,3) = 2*L(1);	D(1,4) = L(2);
	D(2,2) = L(2);		D(2,4) = L(1);		D(2,5) = 2*L(2);
	D(3,3) = 2*L(3);	D(3,4) = L(4);
	D(4,4) = L(3);		D(4,5) = 2*L(4);
	D(5,5) = 2*L(5);
	
	return D;
}



/*******************************************************
****************** Generation Process ******************
********************************************************/

//[[Rcpp::export]]
Matrix_eig Gen_r_from_v_mat(Matrix_eig our_v_mat, Vector_eig sigma){
	int nRow = our_v_mat.rows();	 //n
	int nCol = our_v_mat.cols();	 //m
	Matrix_eig tmp3 = our_v_mat;
	double tmp1, tmp2;
	
	std::srand((unsigned int) time(0));		std::random_device rd{};	std::mt19937 gen{rd()};
	std::normal_distribution<> d{0,1};
	
	for(int i = 0; i < nRow; ++i){		//n
		for(int j = 0; j < nCol; ++j){	//m
			tmp1 = d(rd)*sigma(j);
			tmp2 = d(rd)*sigma(j) + our_v_mat(i,j);		//R::rnorm(0.0, 1.0)
			tmp3(i,j) = std::sqrt(tmp2*tmp2+tmp1*tmp1);
		}
	}
	return(tmp3);
}


//[[Rcpp::export]]
Matrix_eig Gen_r(Matrix_eig W, Vector_eig TE, Vector_eig TR, Vector_eig sigma){
	return(Gen_r_from_v_mat(v_mat(W, TE, TR), sigma));
}


/*******************************************************
********************* Derivatives **********************
********************************************************/


double dee_v_ij_dee_W_ik(Matrix_eig W, Vector_eig TE, Vector_eig TR, int i, int j, int k){
	if(k == 1){
		return( exp(TE(j)*log(W(i, 2))) * (1-exp(TR(j)*log(W(i, 1)))) );
	} else if(k == 2){
		return( -W(i, 0) * TR(j) * exp(TE(j)*log(W(i, 2))) * exp((TR(j)-1)*log(W(i, 1))) );
	} else if(k == 3){
		return( W(i, 0) * TE(j) * exp((TE(j)-1)*log(W(i, 2))) * (1-exp(TR(j)*log(W(i, 1)))) );
	} else {
		return -1000000;
	}
}
// There was a mistake, W(i, 3) would be W(i,2) and so on ... C and R indexing case -- corrected


double simple_dee_v_ij_dee_W_ik(Vector_eig W, Vector_eig TE, Vector_eig TR, int j, int k){
	if(W(1) <= 0){							//Changed
		if(k == 1){
			return( exp(TE(j)*log(W(2))) );
		} else if(k == 2){
			return( -W(0) * TR(j) * exp(TE(j)*log(W(2))) );
		} else if(k == 3){
			return( W(0) * TE(j) * exp((TE(j)-1)*log(W(2))) );
		} else {
			return -10000;
		}
	} else{								// Not doing else if again!
		if(k == 1){
			return( std::exp(TE(j)*std::log(W(2))) * (1-std::exp(TR(j)*std::log(W(1)))) );
		} else if(k == 2){
			return( -W(0) * TR(j) * std::exp(TE(j)*std::log(W(2))) * std::exp((TR(j)-1)*std::log(W(1))) );
		} else if(k == 3){
			return( W(0) * TE(j) * std::exp((TE(j)-1)*std::log(W(2))) * (1-std::exp(TR(j)*std::log(W(1)))) );
		} else {
			return -10000;
		}
	}
}



double dee_2_v_ij_dee_W_ik_dee_W_ik1(Matrix_eig W, Vector_eig TE, Vector_eig TR, int i, int j, int k, int k1){
	if(k == 1 && k1 == 1){
		return 0;
	} else if ((k == 1 && k1 == 2) || (k == 2 && k1 == 1)){
		return ( -TR(j) * exp(TE(j)*log(W(i, 2))) * exp((TR(j)-1)*log(W(i, 1))) );
	} else if ((k == 1 && k1 == 3)||(k == 3 && k1 == 1)){
		return ( TE(j) * exp((TE(j)-1)*log(W(i, 2))) * (1-exp(TR(j)*log(W(i, 1)))) );
	} else if(k == 2 && k1 == 2){
		return( - W(i,0) * TR(j) * (TR(j)-1) * exp(TE(j)*log(W(i,2))) * exp((TR(j)-2)*log(W(i, 1))) );
	} else if((k == 2 && k1 == 3)||(k == 3 && k1 == 2)){
		return( - W(i,0) * TR(j) * TE(j) * exp((TE(j)-1)*log(W(i,1))) * exp((TR(j)-1)*log(W(i, 1))) );
	} else if(k == 3 && k1 == 3){
		return ( W(i,0) * TE(j) * (TE(j)-1) * exp((TE(j)-2)*log(W(i,1))) * (1-exp(TR(j)*log(W(i, 1)))) );
	} else {
		return -1000000;
	}
}

double simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(Vector_eig W, Vector_eig TE, Vector_eig TR, int j, int k, int k1){
	if(k == 1 && k1 == 1){
		return 0;
	} else if ((k == 1 && k1 == 2) || (k == 2 && k1 == 1)){
		return ( -TR(j) * std::exp(TE(j)*std::log(W(2))) * std::exp((TR(j)-1)*std::log(W(1))) );
	} else if ((k == 1 && k1 == 3)||(k == 3 && k1 == 1)){
		return ( TE(j) * std::exp((TE(j)-1)*std::log(W(2))) * (1-std::exp(TR(j)*std::log(W(1)))) );
	} else if(k == 2 && k1 == 2){
		return( - W(0) * TR(j) * (TR(j)-1) * std::exp(TE(j)*std::log(W(2))) * std::exp((TR(j)-2)*std::log(W(1))) );
	} else if((k == 2 && k1 == 3)||(k == 3 && k1 == 2)){
		return( - W(0) * TR(j) * TE(j) * std::exp((TE(j)-1)*std::log(W(1))) * std::exp((TR(j)-1)*std::log(W(1))) );
	} else if(k == 3 && k1 == 3){
		return ( W(0) * TE(j) * (TE(j)-1) * std::exp((TE(j)-2)*std::log(W(1))) * (1-std::exp(TR(j)*std::log(W(1)))) );
	} else {
		return -1000000;
	}
}




Eigen::VectorXf read_sd(char* const sd_file, int our_dim_4){

	// Read variance:
	float output[our_dim_4];
	int n = 0;
	/*
	FILE* fp = fopen(sd_file, "r");
	if (fp == NULL) {
		printf("failed to open file\n");
		exit(EXIT_FAILURE);
	}
	while (fscanf(fp, "%f", &output[n++]) != EOF)
		;
	*/
	Eigen::VectorXf sigma = Eigen::VectorXf::Zero(our_dim_4);
	for(int i = 0; i < our_dim_4; ++i){
		//sigma(i) = 5.; 			//output[i]/20.;
		sigma(i) = 15.;
	}
	Debug0("Read the sd file\n----------------\n");
	
	return sigma;
}


#endif	/* MAIN_HEADER */

