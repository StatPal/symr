/*


Different functions to be used in main files.

Main parts:
Debugging level definitions 
1) Group of Bessel functions
2) Small Misc fn (det, head, check bounds, nan, dim etc)
3) MRF class,
   Defnition of Lambda and related matrices and functions(+Kron), 
   Different Reparametrizations (+Cholesky) 
   Bloch transform
4) Generate a sample matrix when mu and sigma are given
5) Hessian matrix and vector related to Delta method



Changes:

float -> double
Xf -> Xd 
3f -> 3d


Precompile header using

g++ scheme_new.hpp -I /usr/include/eigen3 -O3

- Don't do it now. taking huge gch file. 


*/



//#include <RcppEigen.h>

#ifndef MAIN_HEADER
#define MAIN_HEADER

#include <iostream>
#include <iomanip>
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

typedef std::function<double(const Vector_eig &x)> my_function_type;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::Matrix3d Matrix3d_eig;

#define SQ(x) ((x) * (x))

const int IF_DEBUG = 1;
#define ABS_ERR 1e-10



#define DEBUG_LEVEL_0				//Minimal debugs

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


#define DEBUG_LEVEL_2

#ifdef DEBUG_LEVEL_2
#define Debug2(x) {std::cout << "DEBUG 2: "<< x << "\n";}
#else
#define Debug2(x)
#endif



//#define DEBUG_LEVEL_3

#ifdef DEBUG_LEVEL_3
#define Debug3(x) {std::cout << "DEBUG 3: "<< x << "\n";}
#else
#define Debug3(x)
#endif


#define DEBUG_ANOTHER_LEVEL 0




#define DEBUG_LEVEL_times				// Time debugs only

#ifdef DEBUG_LEVEL_times
#define Debug_time(x) {std::cout << "DEBUG t: "<< x << "\n";}
#else
#define Debug_time(x)
#endif



const Matrix_eig G((Matrix_eig(6,9) << 
  1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1).finished());









/********************************************************
*********************** Eigen ***************************
********************************************************/
/*
* .resize may be useful or reshape
* https://gist.github.com/gocarlos/c91237b02c120c6319612e42fa196d77
*/ 









/********************************************************
******************** Bessel functions *******************
********************************************************/


/*
* Modified Bessel Function of First Kind - Hand-written -  would be slow
* Inputs: x, alpha
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




/* 
* Bessel function: reformatted - need to see what macro can be used
* Inputs: x, nu
*/
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





/*
* Ratio of BesselI(x, 1) and BesselI(x, 0)
* No used now, see later!
*/
double ratio_bessel_10(double x){

	if(x<150){
		return (our_bessel_I(x, 1)/our_bessel_I(x, 0));
	} else {
		return (1.0 - 0.5/x);
	}
}



/*
* Ratio of BesselI(x, 2) and BesselI(x, 0)
*/
double ratio_bessel_20(double x){

	if(x<150){
		return (our_bessel_I(x, 2)/our_bessel_I(x, 0));
		//return (std::cyl_bessel_i(2, x)/std::cyl_bessel_i(0, x));
	} else {
		return (1.0 - 2/x);		// Take one more term
	}
}




/*
* Ratio of BesselI(x, 1) and BesselI(x, 0)
* Copied - give source.
*/
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




/**
* log(BesselI(x, 0))
* Copied - give source.
*/
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
***************** Other small functions *****************
********************************************************/


/* 
* Covariance matrix from a data matrix x
*/
// [[Rcpp::export]]
Matrix_eig Cov_1(Matrix_eig x) {
	int nRows = x.rows();
	Vector_eig colMean = x.colwise().mean();
	// https://eigen.tuxfamily.org/dox/group__TutorialReductionsVisitorsBroadcasting.html
	x.rowwise() -= colMean.transpose();			// or x_cen instead of x
	return x.transpose()*x/(nRows-1);			// or x_cen
}



/*
* Similar to head function in R - helpful for debugging
*/
void show_head(Eigen::MatrixXd W, int n = 10){
	std::cout << "head of the matrix:\n" ;
	for(int i = 0; i < n; ++i){
		std::cout << W.row(i) << "\n";
	}
}


/*
* Similar to head function in R for a vector - helpful for debugging
*/
void show_head_vec(Eigen::VectorXd W, int n = 10, int endLine = 0){
	std::cout << "head of the vector:\t" ;
	for(int i = 0; i < n; ++i){
		std::cout << W(i) << ", ";
	}
	if(endLine) 
		std::cout << "\n";
	else 
		std::cout << "\t";
}



/*
* Mean of rice distribution
*/
double mean_rice(double nu, double sigma){
	double x = - SQ(nu)/(2*SQ(sigma));
	return sigma * std::sqrt(M_PI/2) * std::exp(x/2)*( (1-x)*our_bessel_I(-x/2, 0) - x * our_bessel_I(-x/2, 1)) ;
}



/**
* Change vector to matrix
*/
//[[Rcpp::export]]
Matrix_eig to_matrix(Vector_eig v1, int nrow_in, int ncol_in){
	return(Map<MatrixXd> (v1.data(), nrow_in, ncol_in));
	//https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
	//https://stackoverflow.com/questions/52261389/how-to-convert-an-stdvector-to-a-matrix-in-eigen
	//https://stackoverflow.com/questions/32452739/vector-to-matrix/32475129
}



/**
* Change matrix to vector
*/
Vector_eig to_vector(Matrix_eig v1, int is_transpose=0){
	if(!is_transpose){
		return(Map<VectorXd> (v1.data(), v1.rows()*v1.cols()));
	} else {
		return(Map<VectorXd> (v1.transpose().data(), v1.rows()*v1.cols()));
	}
}




/**
Crude Determinant of a sparse matrix - not needed I guess
*/
// [[Rcpp::export]]
double sp_det_1(SpMat A){
	return MatrixXd(A).determinant();
	//https://stackoverflow.com/questions/15484622/how-to-convert-sparse-matrix-to-dense-matrix-in-eigen/15592295
	//https://stackoverflow.com/questions/13033694/eigen-convert-dense-matrix-to-sparse-one
}


/**
 Vectorized log - needed?
*/
Vector_eig log_vec(Vector_eig x){
	Vector_eig y = x;
	for(int i = 0; i < x.size(); ++i){
		y(i) = std::log(x(i));
	}
	return y;
}

/**
 Needed?
*/
double abs_sum(Vector_eig x){
	double temp = 0.;
	for(int i = 0; i < x.size(); ++i){
		temp += std::fabs(x(i));
	}
	return temp;
}



/**
* Not needed now
*/
// [[Rcpp::export]]
double sp_log_det_2(SpMat B){			// Log determinant
	MatrixXd A = MatrixXd(B);
	SelfAdjointEigenSolver<MatrixXd> es(A);				// Marked was in eigen 2
	return log_vec(es.eigenvalues()).sum();
}



/**
* Not needed now
*/
double sp_log_det_7(SpMat A){			// Log determinant - LU
	Eigen::SparseLU<Eigen::SparseMatrix<double>, COLAMDOrdering<int>> solver;
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



/*
* Log determinant of a matrix
* Not used now - see finding det with Cholesky
*/
// [[Rcpp::export]]
double log_det_2(Matrix_eig B){
	SelfAdjointEigenSolver<MatrixXd> es(B);
	return log_vec(es.eigenvalues()).sum();
}



/*
Log determinant for fixed size(3) matrix
*/
double log_det_3(Matrix3d_eig B){
	SelfAdjointEigenSolver<Matrix3d> es(B);
	return log_vec(es.eigenvalues()).sum();
}




/*
* Log determinant of a 3x3 matrix 
* with Cholesky decomposition to avoid numerical error 
*/
double log_det_3_chol(Matrix3d_eig A){
	//Matrix3d_eig L = to_Cholesky(A);			// Also numerical problems
	Matrix3d_eig L( A.llt().matrixL() );
	// It seems that this does not give error, whereas hand-written algo gives error if L(2,2)^2 ~ 0 sometimes
	
	Vector_eig temp(3);
	temp(0) = L(0,0); temp(1) = L(1,1); temp(2) = L(2,2);
	return 2*log_vec(temp).sum();
}
// Bug - Subrata - forgot the 2 - corrected!



/*
// Because something like this was happening (using to_Cholesky()):
//
//DEBUG 1: Psi_inv:
//   0.115604           0 7.54951e-07
//          0    0.115604    -57.2231
//7.54951e-07    -57.2231       28325
//DEBUG 1: Cholesky of Psi_inv:
//  0.340007          0          0
//         0   0.340007          0
//2.2204e-06     -168.3   0.339462
//DEBUG 1: n*log_det_3(Psi_inv) is nan in Q_star_other_param
//DEBUG 1: Eigenvalues of Psi_inv:
//-0.000799466
//    0.115604
//     28325.1*/



/*
* trace of a sparse matrix
*/
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



/*
* Kroneker product of two dense Matrices
* //https://forum.kde.org/viewtopic.php?f=74&t=50952
*/
Matrix_eig Kron_eig(const Matrix_eig &m1, const Matrix_eig &m2){

	Matrix_eig m3(m1.rows()*m2.rows(), m1.cols()*m2.cols());
	for (int i = 0; i < m1.cols(); i++) {
		for (int j = 0; j < m1.rows(); j++) {
			m3.block(i*m2.rows(), j*m2.cols(), m2.rows(), m2.cols()) =  m1(i,j)*m2;
		}
	}
	return m3;
}

/*
* Kroneker product of two Vectors
*/
Vector_eig Kron_vec_eig(const Vector_eig &m1, const Vector_eig &m2){

	Vector_eig m3(m1.size()*m2.size());
	for (int i = 0; i < m1.size(); i++) {
		m3.segment(i*m2.size(), m2.size()) =  m1(i) * m2;
	}
	return m3;
}


/*
* Kroneker product of two Sparse Matrices
*/
SpMat Kron_Sparse_eig(const SpMat &m1, const SpMat &m2){

	// m1.makeCompressed(); m2.makeCompressed();				// const would not allow makeCompressed
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





/*
* Sparse Identity matrix of size n_x
*/
// [[Rcpp::export]]
SpMat I_n(int n_x){
	SpMat temp(n_x, n_x);				//temp.reserve(n_x);
	temp.setIdentity();
	temp.makeCompressed();
	return(temp);
}



/*
* Sparse J_n matrix of size n_x
* One eigenvalue is 0, hence determinant is 0
*/
// [[Rcpp::export]]
SpMat J_n(int n_x){				// has determinant 0???				// When n_x =1
	SpMat temp(n_x, n_x);
	if(n_x==1){
		temp.coeffRef(0,0) = 1;
	} else if(n_x == 2){
		temp.insert(0,0) = 1;
		temp.insert(0,1) = -1;
		temp.insert(1,0) = -1;
		temp.insert(1,1) = 1;
	} else if(n_x>2){
		temp.reserve(3*n_x-4);
		for(int i = 0; i < n_x-1; i++){
			temp.insert(i, i) = 2.0;
			temp.insert(i,i+1) = -1;
			temp.insert(i+1,i) = -1;
		}
		temp.coeffRef(0,0) = 1;
		temp.insert(n_x-1,n_x-1) = 1;
	}
	temp.makeCompressed();
	return(temp);
}







/*
* Eigen values of J_n
* Be careful of the sorting order.
* https://stackoverflow.com/questions/30188482/sparse-eigenvalues-using-eigen3-sparse
*/
Vector_eig eigenvals_J_n(int n) {
	//Matrix_eig A = Matrix_eig(J_n(n));
	//Eigen::SelfAdjointEigenSolver<Matrix_eig> es(A);
	//return (es.eigenvalues());
	
	
	Vector_eig d_vec = Vector_eig::Zero(n);
	for(int i = 0; i < n; ++i){
		d_vec[i] = 2*(1-std::cos(M_PI*i/n));			// Check, this is oppositely sorted - check
	}
	return d_vec;
}





/*
* MRF_param class:
* Contains matrices and functions which are needed for MRF part calculations.
*/
class MRF_param{

  public:
	int n_x, n_y, n_z, n;
	
	SpMat H_1, H_2, H_3;
	Vector_eig eigenval_1_small, eigenval_2_small, eigenval_3_small;
	Vector_eig one_1, one_2, one_3, one_23, one_12;
	Vector_eig eigenval_1, eigenval_2, eigenval_3, eigens;
	Vector_eig final_vec;
	double deriv_1, deriv_2, deriv_3;
	SpMat Lambda_init;
	SpMat Lambda_init_row;
	Matrix_eig tmp_W_Psi_inv;
	Matrix_eig tmp_Wt_L_W;
	Matrix_eig tmp_Lambda_W;
	Matrix_eig Psi_grad = Matrix_eig::Zero(6, 1);
	
	
	//Constructor:
	// MRF_param(){};
	// Might help in main file MRF - yes it resolved the problem.
	// Without this, the errors were:

/*
scheme_new_OSL_EM_15_GEM.cpp: In instantiation of ‘MRF_optim<T>::MRF_optim(typename cppoptlib::BoundedProblem<T>::TVector) [with T = double; typename cppoptlib::BoundedProblem<T>::TVector = Eigen::Matrix<double, -1, 1>]’:
scheme_new_OSL_EM_15_GEM.cpp:403:48:   required from here
scheme_new_OSL_EM_15_GEM.cpp:222:49: error: no matching function for call to ‘MRF_param::MRF_param()’
  222 |   cppoptlib::BoundedProblem<T>(y_.size()), r2(y_){}
      |                                                 ^
In file included from scheme_new_OSL_EM_15_GEM.cpp:24:
scheme_new_numerical.hpp:684:2: note: candidate: ‘MRF_param::MRF_param(int, int, int)’
  684 |  MRF_param(int n_x, int n_y, int n_z) {
      |  ^~~~~~~~~
scheme_new_numerical.hpp:684:2: note:   candidate expects 3 arguments, 0 provided
scheme_new_numerical.hpp:669:7: note: candidate: ‘MRF_param::MRF_param(const MRF_param&)’
  669 | class MRF_param{
      |       ^~~~~~~~~
scheme_new_numerical.hpp:669:7: note:   candidate expects 1 argument, 0 provided
*/

	// https://stackoverflow.com/questions/18971355/no-matching-function-for-call-to-class  - important
	
	
	
	
	
	MRF_param(int n_x, int n_y, int n_z) {
	
		n = n_x * n_y * n_z;
		H_1 = Kron_Sparse_eig( J_n(n_x), I_n(n_y*n_z));
		H_2 = Kron_Sparse_eig( Kron_Sparse_eig(I_n(n_x), J_n(n_y)), I_n(n_z));
		H_3 = Kron_Sparse_eig( I_n(n_x*n_y), J_n(n_z));
		
		
		eigenval_1_small = eigenvals_J_n(n_x);
		eigenval_2_small = eigenvals_J_n(n_y);
		eigenval_3_small = eigenvals_J_n(n_z);
		
		one_1 = Vector_eig::Ones(n_x);
		one_2 = Vector_eig::Ones(n_y);
		one_3 = Vector_eig::Ones(n_z);
		one_23 = Vector_eig::Ones(n_y*n_z);
		one_12 = Vector_eig::Ones(n_x*n_y);
		
		eigenval_1 = Kron_vec_eig(eigenval_1_small, one_23);
		eigenval_2 = Kron_vec_eig(one_1, Kron_vec_eig(eigenval_2_small, one_3) );
		eigenval_3 = Kron_vec_eig(one_12, eigenval_3_small);
		
		final_vec = Vector_eig::Zero(n);
		eigens    = Vector_eig::Zero(n);
		
		
		Lambda_init = 0.1*H_1 + 0.2*H_2 + 0.3*H_3;
		Lambda_init_row = Lambda_init.row(0);
		tmp_W_Psi_inv = Matrix_eig::Zero(n, 3);
		tmp_Wt_L_W = Matrix_eig::Zero(3, 3);
		tmp_Lambda_W = Matrix_eig::Zero(n, 3);
		
	}
	
	
	
	
	/*
	* Sparse Matrix lambda
	* see also https://stackoverflow.com/questions/38839406/eigen-efficient-kronecker-product or kroneckerProduct
	*/
	SpMat Lambda(const Vector_eig &beta){
		Lambda_init = beta(0)*H_1 + beta(1)*H_2 + beta(2)*H_3;		// Check, is it doable?
		return(Lambda_init);
	}
	
	
	// Specific row of Lambda matrix: Creates some problem.
	/*
	SpMat Lambda_row(Vector_eig beta, int i){
		Lambda_init_row = beta(0)*H_1.row(i) + beta(1)*H_2.row(i) + beta(2)*H_3.row(i);
		return(Lambda_init_row);
	}
	*/
	
	
	Vector_eig MRF_grad_fn(const Matrix_eig &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, int i){
		Vector_eig tmp1(3), tmp2(3), tmp3(3);
		tmp_W_Psi_inv.noalias() = W * Psi_inv;
		tmp1 = H_1.row(i) * tmp_W_Psi_inv;
		tmp2 = H_2.row(i) * tmp_W_Psi_inv;
		tmp3 = H_3.row(i) * tmp_W_Psi_inv;
		Vector_eig MRF_grad = beta(0)*tmp1 + beta(1)*tmp2 + beta(2)*tmp3;
		return MRF_grad;
	}
	
	
	/*
	* log determinant of Lambda(beta, n_x, n_y, n_z)
	*/
	double sp_log_det_specific(const Vector_eig &beta, double thres = 0.000001){
		eigens.noalias() = beta(0) * eigenval_1 + beta(1) * eigenval_2 + beta(2) * eigenval_3;
		double temp = 0.0;
		for(int i = 0; i < eigens.size(); ++i){
			if(eigens(i) > thres){
				temp += std::log(eigens(i));
			}
		}
		return temp;
	}
	
	
	
	
	/*
	* The ratio of eigenvalue sum part of Lambda(beta, n_x, n_y, n_z) for the derivative
	* Depends on the value of k (0, 1, 2 - corresponding to derivative wt beta_1, beta_2, beta_3)
	*/
	double sp_log_inv_specific(const Vector_eig &beta, int k, double thres = 0.000001){
	
		eigens.noalias() =  beta(0) * eigenval_1 + beta(1) * eigenval_2 + beta(2) * eigenval_3;
		
		if(k == 0) {
			final_vec.noalias() = beta(0) * eigenval_1;
		} else if(k == 1) {
			final_vec.noalias() = beta(1) * eigenval_2;
		} else if(k == 2) {
			final_vec.noalias() = beta(2) * eigenval_3;
		}
		for(int i = 1; i < n; ++i){
			final_vec(i) /= eigens(i);
		}
		double temp = (double)final_vec.sum();
		return (temp);
	}
	
	
	
	/* 
	* Numerator of the log likelihood from the MRF part:
	*/
	double MRF_log_likeli_num(const Matrix_eig &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta) {
	
		double tmp2 = -(Psi_inv*W.transpose()*Lambda(beta)*W).trace();
		/*
		int i = 0, j = 0, k = 0;
		double tmp1 = 0.0, tmp2 = 0.0;
		for (k = 0; k < H_1.outerSize(); ++k){
			for (SparseMatrix<double>::InnerIterator it(H_1,k); it; ++it){
				// it.value();
				// it.row();   // row index
				// it.col();   // col index (here it is equal to k)
				i = it.row(); j = it.col();
				tmp1 += it.value() * (W.row(j).transpose() * (Psi_inv * W.row(i)));
			}
		}
		tmp2 += beta(0)*tmp1;
		
		tmp1 = 0.0;
		for (k = 0; k < H_2.outerSize(); ++k){
			for (SparseMatrix<double>::InnerIterator it(H_2,k); it; ++it){
				i = it.row(); j = it.col();
				tmp1 += it.value() * (W.row(j).transpose() * (Psi_inv * W.row(i)));
			}
		}
		tmp2 += beta(1)*tmp1;
		
		tmp1 = 0.0;
		for (k = 0; k < H_3.outerSize(); ++k){
			for (SparseMatrix<double>::InnerIterator it(H_3,k); it; ++it){
				i = it.row(); j = it.col();
				tmp1 += it.value() * (W.row(j).transpose() * (Psi_inv * W.row(i)));
			}
		}
		tmp2 += beta(2)*tmp1;
		*/
		
		return tmp2;
	}
	
	
	
	
	/* 
	* log likelihood from the MRF part.
	*/
	double MRF_log_likeli(const Matrix_eig &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta) {
	
		double tmp_num = MRF_log_likeli_num(W, Psi_inv, beta);
		double likeli_sum = ( tmp_num + 3 * sp_log_det_specific(beta) + 
								n * log_det_3(Psi_inv) - 3 * n * log(2*M_PI) )/2;
		
		return likeli_sum;
	}
	
	
	
	// W'Lambda W matrix:
	Matrix_eig Wt_L_W(const Matrix_eig &W, const Vector_eig &beta){
	
		// SpMat Gamma_inv = Lambda(beta);		// Is it okay to say some sparse mat = some sparse mat? (noalias is not possible)
		Lambda_init = Lambda(beta);
		tmp_Lambda_W.noalias() = Lambda_init*W;	// n x 3
		
		return W.transpose()*tmp_Lambda_W;
	}
	
	
	
	
	
	// Derivative of the likelihood
	Vector_eig MRF_log_likeli_grad(const Matrix_eig &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta) {
	
		// SpMat Gamma_inv = Lambda(beta);		// Is it okay to say some sparse mat = some sparse mat? (noalias is not possible)
		Lambda_init = Lambda(beta);
		
		
		Psi_grad.noalias() = 0.5 * G * to_vector(n * Psi_inv.llt().solve(Matrix3d_eig::Identity(3, 3)) - W.transpose()*Lambda_init*W);
	
		tmp_W_Psi_inv.noalias() = W * Psi_inv; 		// Is this direction faster?
		double beta_x_grad = 1.5*sp_log_inv_specific(beta, 0) - 
		                    	0.5*(W.transpose() * H_1 * tmp_W_Psi_inv).trace();
		double beta_y_grad = 1.5*sp_log_inv_specific(beta, 1) - 
		                    	0.5*(W.transpose() *  H_2 * tmp_W_Psi_inv).trace();
		
		
		Vector_eig grad(8);
		grad.segment(0, 6) = Psi_grad;				// 6x1 matrix I guess
		grad(6) = beta_x_grad; grad(7) = beta_y_grad;
	
		return grad;
	}
	
	
	
	// Destructor:
	~MRF_param(){ }
};











/*
* \nu_{ij} as a mx1 vector from one row of W (and TE, TR)
*/
Eigen::VectorXd Bloch_vec(const Vector_eig &W_row, const Vector_eig &TE, const Vector_eig &TR){

	int m = TE.size();
	Eigen::VectorXd tmp = Eigen::VectorXd::Zero(m);
	for(int j = 0; j < m; ++j) {
		tmp(j) = W_row(0);
	}
	if(W_row(2)>0){
		for(int j = 0; j < m; ++j) {
			tmp(j) *= (1-std::exp(TR(j)*std::log(W_row(1))));
		}
	}
	if(W_row(1)>0){
		for(int j = 0; j < m; ++j) {
			tmp(j) *= (1-std::exp(TR(j)*std::log(W_row(1))));
		}
	}
	/*
	if (W_row(2)>0 & W_row(1)>0) {
		for(int j = 0; j < m; ++j) {
			tmp(j) = W_row(0)*std::exp(TE(j)*std::log(W_row(2)))*
		}
	}
	*/
	return tmp;
}



/*
Write in a good format 
*/
void check_bounds(const Matrix_eig &W, const Vector_eig &lb, const Vector_eig &ub){
	int n = W.rows();
	int count = 0;
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < 3; ++j){
			// Debug2(W(i,j));			// WHY!!!!
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


/*
* Prints dim of a matrix in stdout.
*/
void show_dim(const Matrix_eig &A){
	std::cout << "Dimension of the mat: " << A.rows() << " x " << A.cols() << "\n";
}


/*
* Checks whether any vector(x) is inside proper bounds(lb and ub) or not
*/
void check_bounds_vec(const Vector_eig &x, const Vector_eig &lb, const Vector_eig &ub){

	if(x(0)<lb(0) || x(1)<lb(1) || x(2)<lb(2)){
		//std::cout << "Lower bound crossed initially!";
		if(x(0)<lb(0)) Debug2("1st");
		if(x(1)<lb(1)) Debug2("2nd");
		if(x(2)<lb(2)) Debug2("3rd");
		//bad_bound_1++;
	}
	if(x(0)>ub(0) || x(1)>ub(1) || x(2)>ub(2)){
		//std::cout << "Upper Bound crossed initially!";
		if(x(0)>ub(0)) Debug2("1st");
		if(x(1)>ub(1)) Debug2("2nd");
		if(x(2)>ub(2)) Debug2("3rd");
		//bad_bound_2++;
	}

}


/*
* Checks whether there is NaN or not and prints the location in a matrix
*/
int check_nan(const Matrix_eig &A){
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


/*
* Checks whether there is NaN in W and replaces that element from W_old.
*/
int check_nan_W(Matrix_eig& W, const Matrix_eig& W_old){
	int bad_count = 0;
	int i = 0, j = 0;
	for(i = 0; i < W.rows(); ++i){
		for(j = 0; j < W.cols(); ++j){
			if(std::isnan(W(i, j))){
				Debug0("NAN in location: ("<< i << "," << j<< ") of W!");
				W(i, j) = W_old(i, j);
				bad_count++;
				// return 1;
			}
		}
	}
	return bad_count;
}


/*
* Checks whether there is NaN or not and prints the location in a vector
*/
int check_nan_vec(const Vector_eig &A){
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




/*
* Input: W and TE, TR values
* Output: the whole \nu matrix
*/
// [[Rcpp::export]]
Matrix_eig v_mat(const Matrix_eig &W, const Vector_eig &TE, const Vector_eig &TR){
	int nCol = TE.size();	//m
	int nRow = W.rows();	//n
	Matrix_eig tmp = Matrix_eig::Zero(nRow, nCol);		//check the order
	for(int i = 0; i < nRow; ++i) {
		for(int j = 0; j < nCol; ++j) {
			tmp(i,j) = W(i,0);
			if(W(i, 2)>0){								//changed
				tmp(i,j) *= std::exp(TE(j)*std::log(W(i,2)));
			}
			if(W(i, 1)>0){
				tmp(i,j) *= (1-std::exp(TR(j)*std::log(W(i,1))));
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


/*
* Reparametrization to W from rho, T_1, T_2
* Have not used till now - raw rho, T_1, T_2 aree nowhere used
*/
//[[Rcpp::export]]
Matrix_eig to_W(const Vector_eig &rho, const Vector_eig &T_1, const Vector_eig &T_2){
	Matrix_eig W = Matrix_eig::Zero(rho.size(), 3);
	W.col(0) = rho;				//as<arma::vec>(wrap(rho));
	for(int i = 0; i < rho.size(); ++i){
		W(i, 1) = exp(-1/T_1(i));
		W(i, 2) = exp(-1/T_2(i));
	}
	return(W);
}


/*
* Reparametrize everything to one vector of size 3*n+6+2
* Input: W, Psi_inv(symmetric), beta_x, beta_y
* Output: reparametrized vector (needed for the optimization)
*/
Vector_eig to_param_vec(Matrix_eig W, Matrix3d_eig Psi_inv, double beta_x, double beta_y){
	// to_vector can't handle const
	int n = W.rows();		// check
	Vector_eig temp = Vector_eig::Zero(3*n+6+2);

	//int ind =0;
	//for(int j = 0; j < 3; ++j){
	//	for(int i = 0; i < n; ++i){
	//		temp(ind) = W(i, j);
	//		ind++;
	//	}
	//}

	temp.segment(0,3*n)= to_vector(W); 
	// W.resize(n*3, 1) //vectorise( W, 0 );	//Eigen //http://eigen.tuxfamily.org/dox/group__QuickRefPage.html
	temp(3*n+0) = Psi_inv(0,0); temp(3*n+1) = Psi_inv(0,1); temp(3*n+2) = Psi_inv(0,2);
								temp(3*n+3) = Psi_inv(1,1); temp(3*n+4) = Psi_inv(1,2);
															temp(3*n+5) = Psi_inv(2,2);
	temp(3*n+6) = beta_x;
	temp(3*n+7) = beta_y;

	return temp;
}


/*
* I forgot this - The change is in Psi_inv?
* Possibly it is the reparametrization used for gradient calculation of all param
*/
Vector_eig to_param_vec_grad(const Matrix_eig &W, const Matrix_eig &Psi_inv, double beta_x, double beta_y){
	//const removed
	int n = W.rows();		// check
	Vector_eig temp = Vector_eig::Zero(3*n+6+2);
	temp.segment(0,3*n)= to_vector(W); 
	// W.resize(n*3, 1) //vectorise( W, 0 );	//Eigen //http://eigen.tuxfamily.org/dox/group__QuickRefPage.html
	temp.segment(3*n, 6) = Psi_inv;
	temp(3*n+6) = beta_x;
	temp(3*n+7) = beta_y;

	return temp;
}



/*
* Regaining Symmetric Psi_inv(3x3) from temp_psi(6) vector
*/
Matrix3d_eig to_Psi_inv(const Vector_eig &temp_psi){
	Matrix3d_eig Psi_inv = Matrix3d_eig::Zero(3,3);
	Psi_inv(0,0) = temp_psi(0); Psi_inv(0,1) = temp_psi(1); Psi_inv(0,2) = temp_psi(2);
	Psi_inv(1,0) = temp_psi(1); Psi_inv(1,1) = temp_psi(3); Psi_inv(1,2) = temp_psi(4);
	Psi_inv(2,0) = temp_psi(2); Psi_inv(2,1) = temp_psi(4); Psi_inv(2,2) = temp_psi(5);
	return Psi_inv;
}


/*
* Regaining L matrix from temp_L vector (Cholesky part) -- Wait - this is symmetric - not Cholesky?
*/
/*
Matrix3d_eig to_L_mat(const Vector_eig &temp_L){
	Matrix3d_eig L = Matrix3d_eig::Zero(3,3);
	L(0,0) = temp_L(0); L(0,1) = temp_L(1); L(0,2) = temp_L(2);
	L(1,0) = temp_L(1); L(1,1) = temp_L(3); L(1,2) = temp_L(4);
	L(2,0) = temp_L(2); L(2,1) = temp_L(4); L(2,2) = temp_L(5);
	return L;
}
*/  // Not lower Triangular

Matrix3d_eig to_L_mat(const Vector_eig &temp_L){
	Matrix3d_eig L = Matrix3d_eig::Zero(3,3);
	L(0,0) = temp_L(0); L(0,1) = 0.0;       L(0,2) = 0.0;
	L(1,0) = temp_L(1); L(1,1) = temp_L(3); L(1,2) = 0.0;
	L(2,0) = temp_L(2); L(2,1) = temp_L(4); L(2,2) = temp_L(5);
	return L;
}																// Bug -- corrected


/* 
* Reverse transform - intermediate function: 
* Non-zero elements of Cholesky factor 
* input: L: 3x3 Lower triangular matrix
* Output: 6x1 vector
*/
Vector_eig from_L_mat(const Matrix3d_eig &L) {
	Vector_eig temp_L(6);
	temp_L(0) = L(0,0);
	temp_L(1) = L(1,0);
	temp_L(2) = L(2,0);
	temp_L(3) = L(1,1);
	temp_L(4) = L(2,1);
	temp_L(5) = L(2,2);
	
	return temp_L;
}




/*
* llt (opposite of llt decomposition)
*/
Matrix3d_eig from_Cholesky(const Matrix3d_eig &L){
	return (L*L.transpose());
}



/*
* Not used now - this hand written program creates numerical errors in sqrt part!
* Input : Symmetric 3x3 matrix A
* Output: Cholesky Decomposition of A
* Use Matrix3d_eig L( A.llt().matrixL() );
* See the function log_det_3_chol
*/
Matrix3d_eig to_Cholesky(const Matrix3d_eig &A){
	
	Matrix3d_eig L;
	
	L(0,0) = std::sqrt(A(0,0));		L(0,1) = 0; 								L(0,2) = 0;
	L(1,0) = A(1,0)/L(0,0);			L(1,1) = std::sqrt(A(1,1)-SQ(L(1,0)));		L(1,2) = 0;
	L(2,0) = A(2,0)/L(0,0);			L(2,1) = (A(2,1)-L(2,0)*L(1,0))/L(1,1);		L(2,2) = std::sqrt(A(2,2)-SQ(L(2,0))-SQ(L(2,1)));
	
	return L;
}



/*
* Chain rule for converting gradient with Cholesky parametrization
* Input: 6x1 vector (l = nonzero-vec(L))
* Output: 6x6 Lower Traingular(?) matrix: d LL'/d vec(l) ?
* -- Recheck:
* / 2L_0	L_1		L_2		0		0		0    \
* | 0		L_0		0		2L_1	L_2		0    |
* | 0		0		L_2		0		L_1		2L_2 |			// Mistake in first L_2 -> would be L_0  -- corrected
* | 0		0		0		2L_3	L_4		0    |
* | 0		0		0		0		L_3		2L_4 |
* \ 0		0		0		0		0		2L_5 /
* 
* 
* We have calculated dl_star/d vec_symm(A)
* To caculate: dl_star/d vec_chol(L) = dl_star/d vec_symm(A) * d vec_symm(A)/ d vec_chol(L)
* A = LL'
* vec_symm(A) = [a_00, a_10, a_20, a_11, a_12, a_22]
* vec_chol(L) = [l_00, l_10, l_20, l_11, l_12, l_22]				// changed due to change in symmetry??
* The parameter is: 
*/
Matrix_eig to_grad_Cholesky(const Vector_eig &L){
	
	Matrix_eig D = Matrix_eig::Zero(6, 6);
	
	D(0,0) = 2*L(0);	D(0,1) = L(1);		D(0,2) = L(2);
	D(1,1) = L(0);		D(1,3) = 2*L(1);	D(1,4) = L(2);
	D(2,2) = L(0);		D(2,4) = L(1);		D(2,5) = 2*L(2);		//corrected - Subrata
	D(3,3) = 2*L(3);	D(3,4) = L(4);
	D(4,4) = L(3);		D(4,5) = 2*L(4);
	D(5,5) = 2*L(5);
	
	return D;
}








/*******************************************************
****************** Generation Process ******************
********************************************************/


/*
* Input: \nu matrix and sigma
* Output: Generate a sample r matrix
*/
//[[Rcpp::export]]
Matrix_eig Gen_r_from_v_mat(const Matrix_eig &our_v_mat, const Vector_eig &sigma){
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


/*
* Same function as before with different parametrization
*/
//[[Rcpp::export]]
Matrix_eig Gen_r(const Matrix_eig &W, const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma){
	return(Gen_r_from_v_mat(v_mat(W, TE, TR), sigma));
}










/*******************************************************
********************* Derivatives **********************
********************************************************/


// Not needed now - see next one
double dee_v_ij_dee_W_ik(const Matrix_eig &W, const Vector_eig &TE, const Vector_eig &TR, 
						 int i, int j, int k){
	if(k == 0){
		return( exp(TE(j)*log(W(i, 2))) * (1-exp(TR(j)*log(W(i, 1)))) );
	} else if(k == 1){
		return( -W(i, 0) * TR(j) * exp(TE(j)*log(W(i, 2))) * exp((TR(j)-1)*log(W(i, 1))) );
	} else if(k == 2){
		return( W(i, 0) * TE(j) * exp((TE(j)-1)*log(W(i, 2))) * (1-exp(TR(j)*log(W(i, 1)))) );
	} else {
		return -1000000;
	}
}
// There was a mistake, W(i, 3) would be W(i,2) and so on ... C and R indexing case -- corrected


/**
* d\nu_{ij}/dW_{ik} where i is fixed - j and k are taken as inputs.
*/
double simple_dee_v_ij_dee_W_ik(const Vector_eig &W, const Vector_eig &TE, const Vector_eig &TR, 
								int j, int k){

	if( k != 0 && k != 1 && k != 2){
		Debug0("k is not 0/1/2:" << k);
	}
	double deriv = 0.0;
	
	if(k == 0){
		deriv = 1.0;
	} else if(k == 1){
		deriv = -W(0)*TR(j);
	} else if(k == 2){
		deriv = W(0) * TE(j);
	}
	
	if(W(1) > 0){
		if(k == 0){
			deriv *= (1-std::exp(TR(j)*std::log(W(1))));
		} else if(k == 1){
			deriv *= std::exp((TR(j)-1)*std::log(W(1)));
		} else if(k == 2){
			deriv *= (1-std::exp(TR(j)*std::log(W(1)))) ;
		}
	}
	
	if(W(2) > 0){
		if(k == 0){
			deriv *= std::exp(TE(j)*std::log(W(2)));
		} else if(k == 1){
			deriv *= std::exp(TE(j)*std::log(W(2)));
		} else if(k == 2){
			deriv *= std::exp((TE(j)-1)*std::log(W(2))) ;
		}
	}
	
	return deriv;
	
	/*
	if(k == 1){
		return( std::exp(TE(j)*std::log(W(2))) * (1-std::exp(TR(j)*std::log(W(1)))) );
	} else if(k == 2){
		return( -W(0) * TR(j) * std::exp(TE(j)*std::log(W(2))) * std::exp((TR(j)-1)*std::log(W(1))) );
	} else if(k == 3){
		return( W(0) * TE(j) * std::exp((TE(j)-1)*std::log(W(2))) * (1-std::exp(TR(j)*std::log(W(1)))) );
	} else {
		return -10000;
	}
	*/
}
// There was a problem in indexing - C / R style confusion - corrected 

// Corect the next ones!



// Not needed now - see next one
double dee_2_v_ij_dee_W_ik_dee_W_ik1(const Matrix_eig &W, const Vector_eig &TE, const Vector_eig &TR, 
                                     int i, int j, int k, int k1){

	
	if(k == 0 && k1 == 0){
		return 0;
	} else if ((k == 0 && k1 == 1) || (k == 1 && k1 == 0)){
		return ( -TR(j) * exp(TE(j)*log(W(i, 2))) * exp((TR(j)-1)*log(W(i, 1))) );
	} else if ((k == 0 && k1 == 2)||(k == 2 && k1 == 0)){
		return ( TE(j) * exp((TE(j)-1)*log(W(i, 2))) * (1-exp(TR(j)*log(W(i, 1)))) );
	} else if(k == 1 && k1 == 1){
		return( - W(i,0) * TR(j) * (TR(j)-1) * exp(TE(j)*log(W(i,2))) * exp((TR(j)-2)*log(W(i, 1))) );
	} else if((k == 1 && k1 == 2)||(k == 2 && k1 == 1)){
		return( - W(i,0) * TR(j) * TE(j) * exp((TE(j)-1)*log(W(i,1))) * exp((TR(j)-1)*log(W(i, 1))) );
	} else if(k == 2 && k1 == 2){
		return ( W(i,0) * TE(j) * (TE(j)-1) * exp((TE(j)-2)*log(W(i,1))) * (1-exp(TR(j)*log(W(i, 1)))) );
	} else {
		return -1000000;
	}
}
// There is a problem in indexing - C / R style confusion  - corrected




/**
* d^2\nu_ij/dW_{i, k}dW_{i, k1}
*/
double simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(const Vector_eig &W, const Vector_eig &TE, const Vector_eig &TR, 
                                            int j, int k, int k1){

	if( k != 0 && k != 1 && k != 2){
		Debug0("k is not 0/1/2:" << k);
	}
	if( k1 != 0 && k1 != 1 && k1 != 2){
		Debug0("k1 is not 0/1/2:" << k1);
	}
	
	
	if(k == 0 && k1 == 0){
		return 0;
	} else if ((k == 0 && k1 == 1) || (k == 1 && k1 == 0)){
		return ( -TR(j) * std::exp(TE(j)*std::log(W(2))) * std::exp((TR(j)-1)*std::log(W(1))) );
	} else if ((k == 0 && k1 == 2)||(k == 2 && k1 == 0)){
		return ( TE(j) * std::exp((TE(j)-1)*std::log(W(2))) * (1-std::exp(TR(j)*std::log(W(1)))) );
	} else if(k == 1 && k1 == 1){
		return( - W(0) * TR(j) * (TR(j)-1) * std::exp(TE(j)*std::log(W(2))) * std::exp((TR(j)-2)*std::log(W(1))) );
	} else if((k == 1 && k1 == 2)||(k == 2 && k1 == 1)){
		return( - W(0) * TR(j) * TE(j) * std::exp((TE(j)-1)*std::log(W(1))) * std::exp((TR(j)-1)*std::log(W(1))) );
	} else if(k == 2 && k1 == 2){
		return ( W(0) * TE(j) * (TE(j)-1) * std::exp((TE(j)-2)*std::log(W(1))) * (1-std::exp(TR(j)*std::log(W(1)))) );
	} else {
		return -1000000;
	}
}
// There is a problem in indexing - C / R style confusion -corrected
// ISSUE: The TE and TR are scaled to have derivative just greater than 1
// But here we would need the same - but > 2.










Eigen::VectorXd read_sd(char* const sd_file, int our_dim_4){

	// Read variance:
	double output[our_dim_4];		// warning: ISO C++ forbids variable length array ‘output’ [-Wvla]
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
	Eigen::VectorXd sigma = Eigen::VectorXd::Zero(our_dim_4);
	for(int i = 0; i < our_dim_4; ++i){
		//sigma(i) = 5.; 			//output[i]/20.;
		sigma(i) = 15.;
	}
	Debug0("Read the sd file\n----------------\n");
	
	return sigma;
}





void show_dim_sp(SpMat A){
	std::cout << "Dimension of the mat: " << A.rows() << " x " << A.cols() << "\n";
}








/***************************************************
**************** Information Matrix ****************
****************************************************/



/**
* W is a nx3 matrix.
* Hessian matrix(w.r.t. W) would be 3n x 3n matrix - mostly 0's
* d^2 l* / dW_{i'k'} dW_{ik}
* {i,k} are arranged in a fashion so that k remains like a subdivision under division of i
* i.e., 
* l-th column/row represents: 
	k = l%3; 
	i = (int) l/3;
* and 
	l = 3 * i + k;
* 
* Number of non-zero elements per column: <= 7*3
*/
/*
SpMat Hessian_mat(const Matrix_eig &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                  const Vector_eig &TE, const Vector_eig &TR, 
                  const Vector_eig &sigma, const Matrix_eig &r, 
                  int n_x, int n_y, int n_z){

	
	auto time_1_hess = std::chrono::high_resolution_clock::now();
	Debug1("Hessian calculation started without MRF");
	Matrix_eig v = v_mat(W, TE, TR);
	int n = n_x * n_y * n_z;	//n
	int m = v.cols();			//m
	double temp = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;
	SpMat Gamma_inv = Lambda(beta, n_x, n_y, n_z);
	SpMat W_hess(3*n, 3*n);
	W_hess.reserve( VectorXi::Constant(3*n, 7*3) );
	// Reserve 7*3 non-zero's per column - https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
	Debug1("Hessian matrix allocated");
	
	
	// Diagonal parts //
	
	int i = 0, i1 = 0, k = 0, k1 = 0, j = 0;
	Vector_eig temp_vec(3), temp_vec_1(3), temp_vec_2(3);;
	
	for(i = 0; i < n; ++i) {
		//temp_vec = W.row(i);
		for(k = 0; k < 3; ++k) {
			for(k1 = 0; k1 < 3; ++k1) {
				
				temp = 0.;									// Missed this -- correct this - done!
				for(j = 0; j < m ; ++j) {
					
					tmp2 = r(i,j)/SQ(sigma(j));
					tmp3 = - v(i,j)/SQ(sigma(j)) + tmp2 * besselI1_I0(tmp2 * v(i,j));
					temp += tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1);
					
					tmp2 *= v(i,j);
					tmp3 = (1 + ratio_bessel_20(tmp2) - 2*SQ(besselI1_I0(tmp2)) );
					tmp4 = -1/SQ(sigma(j)) +  0.5*SQ(r(i,j)/SQ(sigma(j)))*tmp3;
					temp += tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
									simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1);
				}
				// W_hess.insert(i+k*n, i+k1*n) = temp;		// old way - not very visually pleasing I guess.
				
				W_hess.insert(3 * i + k, 3 * i + k1) = temp + Gamma_inv.coeff(i, i) * Psi_inv(k, k1);
				//https://stackoverflow.com/questions/42376127/how-to-access-a-specific-row-col-index-in-an-c-eigen-sparse-matrix
			}
		}
	}
	
	
	// Off diagonal(only) parts: //
	
	for(i = 0; i < n; ++i){
		for(i1 = 0; i1 < n; ++i1){
			if(i != i1){
				for(k = 0; k < 3; ++k){
					for(k1 = 0; k1 < 3; ++k1){
						W_hess.insert(3 * i + k, 3 * i1 + k1) = Gamma_inv.coeff(i, i1) * Psi_inv(k, k1);
					}
				}
			}
			
		}
	}
	
	W_hess.makeCompressed();
	
	auto time_2_hess = std::chrono::high_resolution_clock::now();
	auto duration_hess = std::chrono::duration_cast<std::chrono::seconds>(time_2_hess - time_1_hess);
	Debug1("Time taken total loop: " << duration_hess.count() << " seconds\n");
	Debug0("Hessian calculated without MRF");
	
	return(W_hess);	//3nx3n
}
*/




// https://stackoverflow.com/questions/47694725/using-inneriterator
/**
* Calculates the \nu_ij for \nu_ij' \Sigma_ij \nu_ij (i.e., w.r.t. \nu_ij)
* where Sigma is an estimation of variance of (W_ik)_{i,k} 
* 
* i.e., it calculates d\nu_ij/dW_ik
* where j is fixed when we consider just one image and i corresponds to the i-th voxel
* 
* So, it would give a 3n x 1 vector: i.e., d\nu_ij/dW_{i1,k} (confusing notation)
* where the value is 0 when i != i1
* 
* So, we would get a sparse vector(/matrix) with 3 non-zero element of total size 3n x 1
* 
* 0:3 + 3*i - th elements would be non-zero and they are :
* d\nu_ij/dW_{i,0}, d\nu_ij/dW_{i,1}, d\nu_ij/dW_{i,2}
*/
/*
SpVec v_grad(const Matrix_eig &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
             const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, const Matrix_eig &r, 
             int n_x, int n_y, int n_z, int i, int j){


	//SpMat grad(3*n_x*n_y*n_z, 1);
	//grad.reserve(VectorXi::Constant(1, 3));
	
	//for(int i1 = 3*i; i1 < 3 * i + 3; ++i1){
	//	grad.insert(i1, 1) = simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, i1 % 3);
	//}

	SpVec grad(3*n_x*n_y*n_z);
	
	for(int i1 = 3*i; i1 < 3 * i + 3; ++i1){
		grad.insert(i1) = simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, i1 % 3);
	}
	return grad;
}
*/



























	// Export the result to a file:
	//	saveAsBitmap(x, n, argv[1]);
























#endif	/* MAIN_HEADER */

