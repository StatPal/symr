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






To Do:
Increase max iteration in variance estimate



** Numerical fix:
Instead of using besselI, use scaled version from gsl.
Also, in some places, handwritten bessel function was still used which is ineffeicient. 
Those are changed. 

mean_rice function uses that scaled version. 




BUG: If TE is not same as 1:k format,
r[i,] and v_new would not match
i.e., proper index if r is not taken.
-- Sorry - NO BUG. r is not passed. train is passed. Hence completely Okay!

BUG: Lambda * Psi 

Some problem with extra added sigma_j^2 in some later versions, previously it was correct. 

was negative/positive. Hessian had this BUG.

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
#include <fstream>
#include <ostream>



#include <cmath> 	//For bessel fn if cpp17
#include <chrono>


extern "C" {
#include <gsl/gsl_sf_bessel.h>		// Try SCALED besselI from here! - changed a bit
}




// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

//using namespace Rcpp;
using namespace Eigen;


typedef Eigen::MatrixXd Matrix_eig;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, RowMajor> Matrix_eig_row;	// For r, W etc	// Inside optimizer, the THessian would be changed
typedef Eigen::VectorXd Vector_eig;
typedef Eigen::VectorXd::Scalar Scalar_eig;
typedef Eigen::ArrayXd Array_eig;

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


// #define DEBUG_LEVEL_2

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


//#define DEBUG_LEVEL_LS

#ifdef DEBUG_LEVEL_LS
#define DebugLS(x) {std::cout << "DEBUG LS: "<< x << "\n";}
#else
#define DebugLS(x)
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
	
//#ifdef (__cplusplus == 201703L)
	//return(std::cyl_bessel_i(nu, x));		// This was commented most recently
//#else
	//return(besseli(x, nu));
//#endif*/
		
		//return(besseli(x, nu));					//return(R::bessel_i(x, nu, 1));
		
	return(gsl_sf_bessel_In(nu, x));
	
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

	//if(x<150){
		//return (our_bessel_I(x, 1)/our_bessel_I(x, 0));
	//} else {
		return (gsl_sf_bessel_I1_scaled(x)/gsl_sf_bessel_I0_scaled(x));
		//return (1.0 - 0.5/x);
	//}
}



/*
* Ratio of BesselI(x, 2) and BesselI(x, 0)
*/
double ratio_bessel_20(double x){

	//if(x<150){
		//return (our_bessel_I(x, 2)/our_bessel_I(x, 0));
		//return (std::cyl_bessel_i(2, x)/std::cyl_bessel_i(0, x));
	//} else {
		return (gsl_sf_bessel_In_scaled(2, x)/gsl_sf_bessel_I0_scaled(x));
		// return (1.0 - 2/x);		// Take one more term
	//}
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


double var(const Vector_eig &x){
	double tmp = x.sum()/x.size();
	tmp = -SQ(tmp);
	tmp += x.array().square().sum()/x.size();
	return tmp;
}


double abs_dev_mean(const Vector_eig &x){
	double tmp = x.sum()/x.size();
	Vector_eig y = x - Vector_eig::Constant(x.rows(), x.cols(), tmp);
	tmp = y.array().abs().sum()/x.size();
	return tmp;
}



/*
* Similar to head function in R - helpful for debugging
* Used template here.
* Matrix_eig_row or Matrix_eig or similar
* Use show_head<Matrix_eig_row>(W);
*/
/*
template <typename T>
void show_head(const T &W, int n = 10){
	std::cout << "Head of the matrix:\n" ;
	for(int i = 0; i < n; ++i){
		std::cout << W.row(i) << "\n";
	}
}
*/
// https://stackoverflow.com/questions/27687769/use-different-parameter-data-types-in-same-function-c
// https://stackoverflow.com/questions/12331655/function-overloading-vs-function-templates-c -- This explains
// Sorry - Just keep it as overload for now...:
void show_head(const Matrix_eig &W, int n = 10){
	std::cout << "Head of the matrix:\n" ;
	for(int i = 0; i < n; ++i){
		std::cout << W.row(i) << "\n";
	}
}
void show_head(const Matrix_eig_row &W, int n = 10){
	std::cout << "Head of the matrix:\n" ;
	for(int i = 0; i < n; ++i){
		std::cout << W.row(i) << "\n";
	}
}





/*
* Similar to head function in R for a vector - helpful for debugging
*/
void show_head_vec(const Eigen::VectorXd &W, int n = 10, int endLine = 0){
	std::cout << "Head of the vector:\t" ;
	for(int i = 0; i < n; ++i){
		std::cout << W(i) << ", ";
	}
	if(endLine) 
		std::cout << "\n";
	else 
		std::cout << "\t";
}



void show_head_sp(const SpMat &A, int n = 10){
	std::cout << "Head of the Sparse matrix:\n" ;
	for(int i = 0; i < n; ++i){
		std::cout << A.row(i) << "\n";
	}
}




/*
* Mean of rice distribution
* Using GSL for SCALED besselI
*/

double mean_rice(double nu, double sigma){
	double x = - SQ(nu)/(2*SQ(sigma));
	// return sigma * std::sqrt(M_PI/2) * std::exp(x/2)*( (1-x)*our_bessel_I(-x/2, 0) - x * our_bessel_I(-x/2, 1)) ;
	return sigma * std::sqrt(M_PI/2) *( (1-x)*gsl_sf_bessel_I0_scaled(-x/2) - x * gsl_sf_bessel_I1_scaled(-x/2)) ;
	
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
* Used overloading here.
*/
Vector_eig to_vector(Matrix_eig v1, int is_transpose=0){
	if(!is_transpose){
		return(Map<VectorXd> (v1.data(), v1.rows()*v1.cols()));
	} else {
		return(Map<VectorXd> (v1.transpose().data(), v1.rows()*v1.cols()));
	}
}
// Can't be done using Map I guess for rowmajor - ignore transpose now.
// Still has some problem - will see later.
/*
Vector_eig to_vector(Matrix_eig_row v1, int is_transpose=0){
	//if(!is_transpose){
	//	return(Map<VectorXd> (v1.data(), v1.rows()*v1.cols()));
	//} else {
	//	return(Map<VectorXd> (v1.transpose().data(), v1.rows()*v1.cols()));
	//}
	
	// cout << "In memory (row-major):" << endl;
	Vector_eig tmp = Vector_eig::Zero(v1.size());
	for (int i = 0; i < v1.size(); i++){
		tmp(i) = *(v1.data() + i);
	}
	return tmp;
	// or just pass the pointer? Check what is there!
}
*/
Vector_eig to_vector_1(Matrix_eig_row v1, int is_transpose=0){
	//if(!is_transpose){
	//	return(Map<VectorXd> (v1.data(), v1.rows()*v1.cols()));
	//} else {
	//	return(Map<VectorXd> (v1.transpose().data(), v1.rows()*v1.cols()));
	//}
	
	// int n = v1.rows()*v1.cols();
	// Vector_eig tmp(n);
	// Eigen::Map< Vector_eig > tmp(v1.data(), n, 1);
	// https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
	// No idea why it's not working!
	
	Vector_eig tmp = Vector_eig::Zero(v1.size());
	for (int i = 0; i < v1.size(); i++){
		tmp(i) = *(v1.data() + i);
	}
	return tmp;
	// or just pass the pointer? Check what is there!
}
// 2nd part compiled perfectly when name is changed - problem in overloading maybe!
// Though 1st part does not compile even when name is changed!
// https://stackoverflow.com/questions/28722899/creating-an-eigen-matrix-from-an-array-with-row-major-order
// ^ Would this help?




/**
Crude Determinant of a sparse matrix - not needed I guess
*/
// [[Rcpp::export]]
double sp_det_1(const SpMat &A){
	return MatrixXd(A).determinant();
	//https://stackoverflow.com/questions/15484622/how-to-convert-sparse-matrix-to-dense-matrix-in-eigen/15592295
	//https://stackoverflow.com/questions/13033694/eigen-convert-dense-matrix-to-sparse-one
}


/**
 Vectorized log - needed?
*/
Vector_eig log_vec(const Vector_eig &x){
	Vector_eig y = x;
	for(int i = 0; i < x.size(); ++i){
		y(i) = std::log(x(i));
	}
	return y;
}

/**
 Needed?
*/
double abs_sum(const Vector_eig &x){
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
double sp_log_det_2(const SpMat &B){					// Log determinant
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
double log_det_2(const Matrix_eig &B){
	SelfAdjointEigenSolver<MatrixXd> es(B);
	return log_vec(es.eigenvalues()).sum();
}



/*
Log determinant for fixed size(3) matrix
*/
double log_det_3(const Matrix3d_eig &B){
	SelfAdjointEigenSolver<Matrix3d> es(B);
	return log_vec(es.eigenvalues()).sum();
}




/*
* Log determinant of a 3x3 matrix 
* with Cholesky decomposition to avoid numerical error 
*/
double log_det_3_chol(const Matrix3d_eig &A){
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
double sp_trace(const SpMat &A){
	double sum =0;
	for (int k = 0; k < A.outerSize(); ++k) {
		sum += A.coeff(k,k);
	}
	return sum;
}




/*
Write in a good format 
*/
void check_bounds(const Matrix_eig_row &W, const Vector_eig &lb, const Vector_eig &ub){
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
void show_dim(const Matrix_eig_row &A){
	std::cout << "Dimension of the mat: " << A.rows() << " x " << A.cols() << "\n";
}
*/
// It says overload is ambiguous - I guess the second one is not necessary!


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



int check_bounds_vec_3(const Vector_eig &x, const Vector_eig &lb, const Vector_eig &ub){

	int bad_bound_1 = 0, bad_bound_2 = 0;
	if(x(0)<lb(0) || x(1)<lb(1) || x(2)<lb(2)){
		//std::cout << "Lower bound crossed initially!";
		bad_bound_1++;
	}
	if(x(0)>ub(0) || x(1)>ub(1) || x(2)>ub(2)){
		//std::cout << "Upper Bound crossed initially!";
		bad_bound_2++;
	}
	return (bad_bound_1 + bad_bound_2);
}






/*
* Checks whether there is NaN or not and prints the location in a matrix
* returns number of such cases.
*/
int check_nan(const Matrix_eig_row &A, const char* mat_name = ""){
	int bad_count = 0;
	std::cout << mat_name;
	for(int i = 0; i < A.rows(); ++i){
		for(int j = 0; j < A.cols(); ++j){
			if(std::isnan(A(i, j))){
				Debug0("NAN in location: ("<< i << "," << j<< ")!");
				bad_count++;
			}
		}
	}
	return bad_count;
}


/*
* Checks whether there is NaN in W and replaces that element from W_old.
* Subrata - Important - Replace W row wise - i.e., when replace, replace the whole row
* OW there is no gurentee that it does not at least decrease the log likelihood!
*/
int check_nan_W(Matrix_eig_row& W, const Matrix_eig_row& W_old){
	int bad_count = 0;
	int i = 0, j = 0;
	for(i = 0; i < W.rows(); ++i){
		for(j = 0; j < W.cols(); ++j){
			if(std::isnan(W(i, j))){
				Debug0("NAN in location: ("<< i << "," << j<< ") of W!");
				// W(i, j) = W_old(i, j);		// Subrata - see what to do...
				W.row(i) = W_old.row(i);		// Is this better?
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
	int bad_count = 0;
	for(int i = 0; i < A.size(); ++i){
		if(std::isnan(A(i))){
			Debug0("NAN in location: ("<< i << ") of the vector!" << std::flush);
			bad_count++;
		}
	}
	return bad_count;
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
* Be careful of the sorting order. This is oppositely sorted
* https://stackoverflow.com/questions/30188482/sparse-eigenvalues-using-eigen3-sparse
*/
Vector_eig eigenvals_J_n(int n) {
	//Matrix_eig A = Matrix_eig(J_n(n));
	//Eigen::SelfAdjointEigenSolver<Matrix_eig> es(A);
	//return (es.eigenvalues());
	
	Vector_eig d_vec = Vector_eig::Zero(n);
	for(int i = 0; i < n; ++i){
		d_vec[i] = 2*(1-std::cos(M_PI*i/n));
	}
	return d_vec;
}





/*
* MRF_param class:
* Contains matrices and functions which are needed for MRF part calculations.
*/
class MRF_param{

  public:
	int n_x_, n_y_, n_z_, n;
	
	SpMat H_1, H_2, H_3;
	Vector_eig eigenval_1_small, eigenval_2_small, eigenval_3_small;
	Vector_eig one_1, one_2, one_3, one_23, one_12;
	Vector_eig eigenval_1, eigenval_2, eigenval_3, eigens;
	Vector_eig final_vec;
	double deriv_1, deriv_2, deriv_3;
	SpMat Lambda_init;
	SpMat Lambda_init_row;
	SpMat tmp_row;
	Matrix_eig tmp_W_Psi_inv;
	Matrix_eig tmp_Wt_L_W;
	Matrix_eig tmp_Lambda_W;
	Matrix_eig Psi_grad = Matrix_eig::Zero(6, 1);
	// Matrix_eig tmp_i_Psi_inv = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_Psi_inv_new = Matrix_eig::Zero(1, 3);
	
	
	
	// Important variables:
	Vector_eig x_old;
	double old_MRF_likeli = 0.0;
	double tmp_i_coeff_1 = 0.0;
	
	
	
	// For likeli num i
	Matrix_eig tmp_i_1 = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_2 = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_3 = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_Psi_inv_1 = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_Psi_inv_2 = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_Psi_inv_3 = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_Psi_inv_final = Matrix_eig::Zero(1, 3);
	
	
	Vector_eig tmp1_vec = Vector_eig::Zero(3);
	Vector_eig tmp2_vec = Vector_eig::Zero(3);
	Vector_eig tmp3_vec = Vector_eig::Zero(3);
	Vector_eig MRF_grad = Vector_eig::Zero(3);
	
	
	
	//Constructor:
	// MRF_param(){};
	// https://stackoverflow.com/questions/18971355/no-matching-function-for-call-to-class  - important
	
	
	
	
	
	MRF_param(int n_x, int n_y, int n_z) {
	
		n_x_ = n_x; n_y_ = n_y; n_z_ = n_z;
		n = n_x * n_y * n_z;
		/*
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
		*/
		
		
		H_1 = Kron_Sparse_eig( I_n(n_z*n_y), J_n(n_x));
		H_2 = Kron_Sparse_eig( Kron_Sparse_eig(I_n(n_z), J_n(n_y)), I_n(n_x));
		H_3 = Kron_Sparse_eig( J_n(n_z), I_n(n_y*n_x));
		
		eigenval_1_small = eigenvals_J_n(n_x);
		eigenval_2_small = eigenvals_J_n(n_y);
		eigenval_3_small = eigenvals_J_n(n_z);
		
		one_1 = Vector_eig::Ones(n_x);
		one_2 = Vector_eig::Ones(n_y);
		one_3 = Vector_eig::Ones(n_z);
		one_23 = Vector_eig::Ones(n_y*n_z);
		one_12 = Vector_eig::Ones(n_x*n_y);
		
		eigenval_1 = Kron_vec_eig(one_23, eigenval_1_small);
		eigenval_2 = Kron_vec_eig(one_3, Kron_vec_eig(eigenval_2_small, one_1) );
		eigenval_3 = Kron_vec_eig(eigenval_3_small, one_12);
		
		
		
		
		
		final_vec = Vector_eig::Zero(n);
		eigens    = Vector_eig::Zero(n);
		
		
		Lambda_init = 0.1*H_1 + 0.2*H_2 + 0.3*H_3;
		Lambda_init_row = Lambda_init.row(0);
		tmp_row = H_1.row(0);
		tmp_W_Psi_inv = Matrix_eig::Zero(n, 3);
		tmp_Wt_L_W = Matrix_eig::Zero(3, 3);
		tmp_Lambda_W = Matrix_eig::Zero(n, 3);
		x_old = Vector_eig::Zero(3);
		MRF_grad = Vector_eig::Zero(3);
	}
	
	
	
	
	/*
	* Sparse Matrix lambda
	* see also https://stackoverflow.com/questions/38839406/eigen-efficient-kronecker-product or kroneckerProduct
	*/
	SpMat Lambda(const Vector_eig &beta){
		Lambda_init = beta(0)*H_1 + beta(1)*H_2 + beta(2)*H_3;		// Check, is it doable?
		return(Lambda_init);
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
			final_vec.noalias() = eigenval_1;
			// final_vec.noalias() = beta(0) * eigenval_1;
		} else if(k == 1) {
			final_vec.noalias() = eigenval_2;
			// final_vec.noalias() = beta(1) * eigenval_2;
		} else if(k == 2) {
			final_vec.noalias() = eigenval_3;
			// final_vec.noalias() = beta(2) * eigenval_3;
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
	double MRF_log_likeli_num(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta) {
	
		//double tmp2 = (Psi_inv * W.transpose() * Lambda(beta) * W).trace();
		
		int k = 0;
		double tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;
		for (k = 0; k < H_1.outerSize(); ++k){
			for (SpMat::InnerIterator it(H_1,k); it; ++it){
				// it.value();
				// i = it.row();   // row index
				// j = it.col();   // col index (here it is equal to k)
				tmp3 = (W.row(it.col()) * (Psi_inv*W.row(it.row()).transpose()));
				tmp1 += tmp3 * it.value();
			}
		}
		tmp2 += beta(0)*tmp1;
		
		tmp1 = 0.0;
		for (k = 0; k < H_2.outerSize(); ++k){
			for (SpMat::InnerIterator it(H_2,k); it; ++it){
				tmp3 = (W.row(it.col()) * (Psi_inv*W.row(it.row()).transpose()));
				tmp1 += tmp3 * it.value();
			}
		}
		tmp2 += beta(1)*tmp1;
		
		tmp1 = 0.0;
		for (k = 0; k < H_3.outerSize(); ++k){
			for (SpMat::InnerIterator it(H_3,k); it; ++it){
				tmp3 = (W.row(it.col()) * (Psi_inv*W.row(it.row()).transpose()));
				tmp1 += tmp3 * it.value();
			}
		}
		tmp2 += beta(2)*tmp1;
		
		// Debug0("tmp2: " << tmp2);
		return (-tmp2/2);		// -0.5*trace(Psi_inv*W.transpose()*Lambda(beta)*W);
	}
	
	
	
	
	/* 
	* log likelihood from the MRF part.
	*/
	double MRF_log_likeli(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta) {
	
		double likeli_sum = MRF_log_likeli_num(W, Psi_inv, beta);
		// Debug0("likeli_sum num: " << likeli_sum);
		likeli_sum += ( 3 * sp_log_det_specific(beta) + 
								n * log_det_3(Psi_inv) - 3 * n * log(2*M_PI) )/2;
		//Debug0("likeli_sum den: "<< ( 3 * sp_log_det_specific(beta) + 
		//						n * log_det_3(Psi_inv) - 3 * n * log(2*M_PI) )/2  );
		
		return likeli_sum;
		
		// (-0.5*trace(Psi_inv*W.transpose()*Lambda(beta)*W)) 
		// 			- 1.5*n*log(2*M_PI) + 1.5*log|Gamma^{-1}| + (n/2)*log|Psi^{-1}|
	}
	
	
	
	
	
	
	
	/*
	* Values are updated: 
	* To find 
			sum_{j != i} (Lambda(i, j) * W.row(j)) * Psi_inv 				// Not exactly this -- BUG
	*
	* splitted in two parts:
			sum_{j != i} (H_1(i, j) * W.row(j)) * Psi_inv * beta(0)			// BUG this part * 2 
		and
			sum_{j != i} (H_2(i, j) * W.row(j)) * Psi_inv 
	*/
	void update_neighbours_likeli_old(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, 
								      const Vector_eig &beta, const int i){
	
		tmp_i_1 = Matrix_eig::Zero(1, 3);
		tmp_i_2 = Matrix_eig::Zero(1, 3);
		tmp_i_3 = Matrix_eig::Zero(1, 3);
		
		tmp_row = H_1.row(i);							//Check this can be omitted or not - might take some time
		for(int k = 0; k < tmp_row.outerSize(); ++k){
			for (SpMat::InnerIterator it(tmp_row, k); it; ++it){
				if(it.col() != i){
					tmp_i_1 += W.row(it.col()) * it.value();		//* beta(0)
				}
			}
		}
		tmp_i_1 = tmp_i_1 * beta(0);
		
		
		//if(i == 99){
			//Debug0("tmp_i_1: " << tmp_i_1);
			//Debug0("Check 0.1.1, i = " << i );
		//}
		
		tmp_row = H_2.row(i);
		for(int k = 0; k < tmp_row.outerSize(); ++k){
			for (SpMat::InnerIterator it(tmp_row, k); it; ++it){
				if(it.col() != i){
					tmp_i_2 += W.row(it.col()) * it.value();
				}
			}
		}
		tmp_i_2 = tmp_i_2 * beta(1);
		
		
		
		tmp_row = H_3.row(i);
		for(int k = 0; k < tmp_row.outerSize(); ++k){
			for (SpMat::InnerIterator it(tmp_row, k); it; ++it){
				if(it.col() != i){
					tmp_i_3 += W.row(it.col()) * it.value();
				}
			}
		}
		
		/*
		tmp_i_Psi_inv_1.noalias() = tmp_i_1 * Psi_inv;
		tmp_i_Psi_inv_2.noalias() = tmp_i_2 * Psi_inv;
		tmp_i_Psi_inv_3.noalias() = tmp_i_3 * Psi_inv;
		tmp_i_Psi_inv_final.noalias() = tmp_i_Psi_inv_1 + tmp_i_Psi_inv_2 + tmp_i_Psi_inv_3; //Was not multiplied by 2 -- BUG
		tmp_i_Psi_inv_final *= 2;													// BUG fixed I guess
		*/
		
		tmp_i_Psi_inv_final.noalias() = tmp_i_1 + tmp_i_2 + tmp_i_3;
		tmp_i_Psi_inv_final = tmp_i_Psi_inv_final * Psi_inv;
		tmp_i_Psi_inv_final *= 2;													// BUG fixed I guess
		
		tmp_i_coeff_1 = beta(0) * H_1.coeff(i, i) + beta(1) * H_2.coeff(i, i) + H_3.coeff(i, i);
	}
	
	
	
	
	
	
	/*
	* Values are updated: 
	* To find 
			sum_{j != i} (Lambda(i, j) * W.row(j)) * Psi_inv 				// Not exactly this -- BUG
	*
	* splitted in two parts:
			sum_{j != i} (H_1(i, j) * W.row(j)) * Psi_inv * beta(0)			// BUG this part * 2 
		and
			sum_{j != i} (H_2(i, j) * W.row(j)) * Psi_inv 
	*/
	void update_neighbours_likeli(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, 
								  const Vector_eig &beta, const int i){
	
		tmp_i_1 = Matrix_eig::Zero(1, 3);
		tmp_i_2 = Matrix_eig::Zero(1, 3);
		tmp_i_3 = Matrix_eig::Zero(1, 3);
		
		for (SpMat::InnerIterator it(H_1, i); it; ++it){
			if(it.row() != i){
				tmp_i_1 += W.row(it.row()) * it.value();		//* beta(0)
			}
		}
		tmp_i_1 = tmp_i_1 * beta(0);
		
		
		
		for (SpMat::InnerIterator it(H_2, i); it; ++it){
			if(it.row() != i){
				tmp_i_2 += W.row(it.row()) * it.value();
			}
		}
		tmp_i_2 = tmp_i_2 * beta(1);
		
		
		
		for (SpMat::InnerIterator it(H_3, i); it; ++it){
			if(it.row() != i){
				tmp_i_3 += W.row(it.row()) * it.value();
			}
		}
		
		
		/*
		tmp_i_Psi_inv_1.noalias() = tmp_i_1 * Psi_inv;
		tmp_i_Psi_inv_2.noalias() = tmp_i_2 * Psi_inv;
		tmp_i_Psi_inv_3.noalias() = tmp_i_3 * Psi_inv;
		tmp_i_Psi_inv_final.noalias() = tmp_i_Psi_inv_1 + tmp_i_Psi_inv_2 + tmp_i_Psi_inv_3; //Was not multiplied by 2 -- BUG
		tmp_i_Psi_inv_final *= 2;													// BUG fixed I guess
		*/
		
		tmp_i_Psi_inv_final.noalias() = tmp_i_1 + tmp_i_2 + tmp_i_3;
		tmp_i_Psi_inv_final = tmp_i_Psi_inv_final * Psi_inv;
		tmp_i_Psi_inv_final *= 2;													// BUG fixed I guess
		
		tmp_i_coeff_1 = beta(0) * H_1.coeff(i, i) + beta(1) * H_2.coeff(i, i) + H_3.coeff(i, i);
	}
	
	
	
	
	
	
	
	
	
	
	/* 
	* Numerator of the log likelihood from the MRF part:
	* w.r.t. i-th row of W.
		tmp_i_Psi_inv_new = sum_j (Lambda(i, j) * W.row(j)) * Psi_inv 					// Not this -- BUG
								 = sum_{j != i} (Lambda(i, j) * W.row(j)) * Psi_inv + 	// This line multiplied by 2
								 + Lambda(i, i) * x' * Psi_inv
	* final output is: 
		 - tmp_i_Psi_inv_new * x / 2;
	*/
	double MRF_log_likeli_num_i_new(const Vector_eig &x, const Matrix3d_eig &Psi_inv) {
	
		//if(i == 99){
		//	Debug0("W.row(i): " << W.row(i));
		//	Debug0("(beta(0) * H_1.coeff(i, i) ) * x.transpose(): " << (beta(0) * H_1.coeff(i, i) ) * x.transpose());
		//	Debug0("new MRF log likeli i:" << -0.5 * (tmp_i_Psi_inv_new * x).value());
		//}
		
		//tmp_i_Psi_inv_new.noalias() = tmp_i_Psi_inv_final + 
		//								tmp_i_coeff_1 * (x.transpose() * Psi_inv);		// Can be little more compact here
		tmp_i_Psi_inv_new.noalias() = (x.transpose() * Psi_inv);
		tmp_i_Psi_inv_new = tmp_i_Psi_inv_final + tmp_i_coeff_1 * tmp_i_Psi_inv_new;
		
		return ( -0.5 * (tmp_i_Psi_inv_new * x).value());
	}
	
	
	
	
	
	
	
	
	
	
	
	/*
	* gradient of likelihood w.r.t. i-th row of W.
	* Shorten if possible - maybe using increment
	*/
	Vector_eig MRF_grad_fn_old(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, int i){
		
		tmp_W_Psi_inv.noalias() = W * Psi_inv;
		tmp1_vec.noalias() = H_1.row(i) * tmp_W_Psi_inv;
		tmp2_vec.noalias() = H_2.row(i) * tmp_W_Psi_inv;
		tmp3_vec.noalias() = H_3.row(i) * tmp_W_Psi_inv;
		MRF_grad.noalias() = beta(0)*tmp1_vec + beta(1)*tmp2_vec + beta(2)*tmp3_vec;
		return MRF_grad;
	}
	
	
	
	
	void MRF_grad_fn(const Vector_eig &x, const Matrix3d_eig &Psi_inv, Vector_eig &MRF_grad){
	
		MRF_grad = Psi_inv * x;
		MRF_grad *= tmp_i_coeff_1;
		MRF_grad += 0.5 * tmp_i_Psi_inv_final.transpose();
	}

	
	
	
	
	
	/* 
	* Change in Numerator of the log likelihood from the MRF part
	* if one changes i-th row of W.
	* l(W_new) - l(W_old) 		// In the notebook, -ve of this is calculated.
	* Directly using x_old is better than W_old because W_old might be already changed in the previous loop.
	*/
	double MRF_log_likeli_num_through_increment(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
												int i1) {
	
		Vector_eig x_new = W.row(i1);
		
		// int i = 0, j = 0, k = 0;
		double tmp1 = 0.0, tmp3 = 0.0;
		// double tmp2 = 0.0;
		Lambda_init = Lambda(beta);								// Remove later and split into H_i's
		
		
		// here k and it.col() are same;
		// So we need 
		// \Sum_{j \ne i1} Lambda(i1, j) [W_j * Psi_inv * ( W_i1_new - W_i1_old)']
		// 		+ 0.5 Lambda(i1, i1) [W_i1_new * Psi_inv * W_i1_new' - W_i1_old * Psi_inv * W_i1_old']
		
		// it.col() is fixed at i1.
		// it.row() is basically j
		for (SpMat::InnerIterator it(Lambda_init, i1); it; ++it){
			if(it.col() != i1){
				tmp3 = (W.row(it.row()) * Psi_inv) * (W.row(i1) - x_old).transpose();
				tmp1 += tmp3 * it.value();
			} else if(it.col() == i1){
				tmp3 = (W.row(i1) - x_old) * Psi_inv * (W.row(i1) - x_old).transpose();
				tmp1 += 0.5 * tmp3 * it.value();
			}
		}
		tmp1 = -tmp1;	// = l(W_new) - l(W_old)	// In the notebook, -ve of this is calculated.
		
		return (old_MRF_likeli + tmp1);
	}
	
	
	
	
	
	
	
	// W'Lambda W matrix:
	Matrix_eig Wt_L_W(const Matrix_eig_row &W, const Vector_eig &beta){
	
		// SpMat Gamma_inv = Lambda(beta);
		// Is it okay to say some sparse mat = some sparse mat? (noalias is not possible)
		Lambda_init = Lambda(beta);
		tmp_Lambda_W.noalias() = Lambda_init*W;
		
		return W.transpose()*tmp_Lambda_W;
	}
	
	
	
	
	
	// Derivative of the MRF-likelihood w.r.t. MRF parameters
	Vector_eig MRF_log_likeli_grad(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta) {
	
		// SpMat Gamma_inv = Lambda(beta);
		// Is it okay to say some sparse mat = some sparse mat? (noalias is not possible)
		Lambda_init = Lambda(beta);
		
		
		Psi_grad.noalias() = 0.5 * G * to_vector(n * Psi_inv.llt().solve(Matrix3d_eig::Identity(3, 3)) - W.transpose()*Lambda_init*W);
	
		tmp_W_Psi_inv.noalias() = W * Psi_inv; 		// Is this direction faster?, make another temp mat
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
* \nu_{ij} = W_i0 * (1 - W_i1 ^ TR_j ) * W_i2 ^TE_j
*/
/*
Eigen::VectorXd Bloch_vec(const Vector_eig &W_row, const Vector_eig &TE, const Vector_eig &TR){

	int m = TE.size();
	Eigen::VectorXd tmp(m);
	for(int j = 0; j < m; ++j) {
		tmp(j) = W_row(0) * 
					(1 - std::exp(TR(j)*std::log(W_row(1)))) *
					std::exp(TE(j)*std::log(W_row(2)));
	}
	return tmp;
}
*/
// pass the vector 
void Bloch_vec(const Vector_eig &W_row, const Vector_eig &TE, const Vector_eig &TR, Vector_eig &tmp){

	int m = TE.size();
	for(int j = 0; j < m; ++j) {
		tmp(j) = W_row(0) * 
					(1 - std::exp(TR(j)*std::log(W_row(1)))) *
					std::exp(TE(j)*std::log(W_row(2)));
	}
}





/*
* Input: W and TE, TR values
* Output: the whole \nu matrix
*/
// [[Rcpp::export]]
Matrix_eig_row v_mat(const Matrix_eig_row &W, const Vector_eig &TE, const Vector_eig &TR){
	int nCol = TE.size();	//m
	int nRow = W.rows();	//n
	Matrix_eig_row tmp = Matrix_eig_row::Zero(nRow, nCol);		//check the order
	
	// Keep an eye: there is another pragma outside in the optimizer
	// I guess that would not be a problem
	#pragma omp parallel for
	for(int i = 0; i < nRow; ++i) {
		for(int j = 0; j < nCol; ++j) {
			tmp(i,j) = W(i,0) * 
						(1-std::exp(TR(j)*std::log(W(i,1)))) * 
						std::exp(TE(j)*std::log(W(i,2)));
			
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
* Have not used till now - raw rho, T_1, T_2 are nowhere used
* SHOULD BE CHANGED AS TE_SCALE AND TR_SCALE HAS ARRIVED.
*/
//[[Rcpp::export]]
Matrix_eig_row to_W(const Vector_eig &rho, const Vector_eig &T_1, const Vector_eig &T_2){
	Matrix_eig_row W = Matrix_eig_row::Zero(rho.size(), 3);
	W.col(0) = rho;
	for(int i = 0; i < rho.size(); ++i){
		W(i, 1) = exp(-1/T_1(i));
		W(i, 2) = exp(-1/T_2(i));
	}
	return(W);
}


/*
* Reparametrize everything to one vector of size 3*n+6+2
* Input: W, Psi_inv(symmetric), beta_x, beta_y
* Output: reparametrized vector (needed for the full set optimization)
* Serious change is needed.
*/
Vector_eig to_param_vec(Matrix_eig W, Matrix3d_eig Psi_inv, double beta_x, double beta_y){
	// to_vector can't handle const - use const and hard code temp -- Subrata 
	// Or, just use W.data() -- I guess it would work - Check
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
// Change this for Matrix_eig_row.





/*
* I forgot this - The change is in Psi_inv?
* Possibly it is the reparametrization used for gradient calculation of all param
*/
Vector_eig to_param_vec_grad(const Matrix_eig_row &W, const Matrix_eig &Psi_inv, double beta_x, double beta_y){
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
* Regaining (lower Triangular) L matrix from temp_L vector (Cholesky part)
*/
Matrix3d_eig to_L_mat(const Vector_eig &temp_L){
	Matrix3d_eig L = Matrix3d_eig::Zero(3,3);
	L(0,0) = temp_L(0); L(0,1) = 0.0;       L(0,2) = 0.0;
	L(1,0) = temp_L(1); L(1,1) = temp_L(3); L(1,2) = 0.0;
	L(2,0) = temp_L(2); L(2,1) = temp_L(4); L(2,2) = temp_L(5);
	return L;
}



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
* So default Cholesky decomposition function is used.
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
* | 0		0		L_0		0		L_1		2L_2 |
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
	D(2,2) = L(0);		D(2,4) = L(1);		D(2,5) = 2*L(2);
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
Matrix_eig_row Gen_r_from_v_mat(const Matrix_eig_row &our_v_mat, const Vector_eig &sigma){
	int nRow = our_v_mat.rows();	 //n
	int nCol = our_v_mat.cols();	 //m
	Matrix_eig_row tmp3 = our_v_mat;
	double tmp1, tmp2;
	
	std::srand((unsigned int) time(0));		std::random_device rd{};	std::mt19937 gen{rd()};
	std::normal_distribution<> d{0,1};
	
	for(int i = 0; i < nRow; ++i){		//n
		for(int j = 0; j < nCol; ++j){	//m
			tmp1 = d(rd)*sigma(j);
			tmp2 = d(rd)*sigma(j) + our_v_mat(i,j);		//R::rnorm(0.0, 1.0)
			tmp3(i,j) = std::sqrt(SQ(tmp2)+SQ(tmp1));
		}
	}
	return(tmp3);
}


/*
* Same function as before with different parametrization
*/
//[[Rcpp::export]]
Matrix_eig_row Gen_r(const Matrix_eig_row &W, const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma){
	return(Gen_r_from_v_mat(v_mat(W, TE, TR), sigma));
}










/*******************************************************
********************* Derivatives **********************
********************************************************/


// Not needed now - see next one
double dee_v_ij_dee_W_ik(const Matrix_eig_row &W, const Vector_eig &TE, const Vector_eig &TR, 
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
		deriv = 1.0 * 
				std::exp(TE(j)*std::log(W(2))) * 
				(1 - std::exp(TR(j)*std::log(W(1))));
	} else if(k == 1){
		deriv = -W(0)*TR(j) * 
				std::exp(TE(j)*std::log(W(2))) * 
				std::exp((TR(j)-1)*std::log(W(1)));
	} else if(k == 2){
		deriv = W(0) * TE(j) *
				std::exp((TE(j)-1)*std::log(W(2))) * 
				(1 - std::exp(TR(j)*std::log(W(1))));
	}
	

	return deriv;
}
// Put back the modifications done to avoid numerical errors.





// Not needed now - see next one
double dee_2_v_ij_dee_W_ik_dee_W_ik1(const Matrix_eig_row &W, const Vector_eig &TE, const Vector_eig &TR, 
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
		return( - W(i,0) * TR(j) * TE(j) * exp((TE(j)-1)*log(W(i,2))) * exp((TR(j)-1)*log(W(i, 1))) );
	} else if(k == 2 && k1 == 2){
		return ( W(i,0) * TE(j) * (TE(j)-1) * exp((TE(j)-2)*log(W(i,2))) * (1-exp(TR(j)*log(W(i, 1)))) );
	} else {
		return -1000000;
	}
}
// BUG, see next




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
		return( - W(0) * TR(j) * TE(j) * std::exp((TE(j)-1)*std::log(W(2))) * std::exp((TR(j)-1)*std::log(W(1))) );
	} else if(k == 2 && k1 == 2){
		return ( W(0) * TE(j) * (TE(j)-1) * std::exp((TE(j)-2)*std::log(W(2))) * (1-std::exp(TR(j)*std::log(W(1)))) );
	} else {
		return -1000000;
	}
}
// BUG in last two cases.
// W(1) instead of W(1) 
// ISSUE: The TE and TR are scaled to have derivative just greater than 1
// But here we would need the same - but > 2.










Eigen::VectorXd read_sd(char* const sd_file, int our_dim_4){

	// Read sd:
	// double output[our_dim_4];		// warning: ISO C++ forbids variable length array ‘output’ [-Wvla]
	// no need 							// Change this
	
	
	
	
	std::fstream myfile(sd_file, std::ios_base::in);
	Eigen::VectorXd sigma = Eigen::VectorXd::Zero(our_dim_4);

	float a;
	int  i = 0;
	while (myfile >> a) {
		sigma(i) = a;
		i++;
	}
	
	/*
	FILE* fp = fopen(sd_file, "r");
	if (fp == NULL) {
		printf("failed to open file\n");
		exit(EXIT_FAILURE);
	}
	while (fscanf(fp, "%f", &output[n++]) != EOF)
		;
	*/
	/*
	for(int i = 0; i < our_dim_4; ++i){
		//sigma(i) = 5.; 			//output[i]/20.;
		sigma(i) = 15.;
	}
	*/
	Debug0("Read the sd file\n----------------\n");
	
	return sigma;
}





void show_dim_sp(SpMat A){
	std::cout << "Dimension of the mat: " << A.rows() << " x " << A.cols() << "\n";
}








/***************************************************
**************** Information Matrix ****************
****************************************************/


/*
* Derivative of I_1/I_0
*/
double h(double x){
	double tmp = (1.0 + ratio_bessel_20(x) - 2*SQ(ratio_bessel_10(x)) ); // besselI1_I0 replaced
	return(0.5*tmp);
}








// Compressed Column (or Row) Storage schemes (CCS or CRS)
// https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
// https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html

void save_sparse(Eigen::SparseMatrix<double> sm, const char* file_name, int if_triplet
				// , int if_nnz = 1
				){
	
	
	// _Scalar is double here
	// _StorageIndex default is int
	// change if necessary.
	
	
	std::ofstream file_connection;
	file_connection.open(file_name);
	// Use std::setprecision(8) if necessary

	sm.makeCompressed();
	
	
	
	if(if_triplet){
		for (int k = 0; k < sm.outerSize(); ++k) {
			for (SpMat::InnerIterator it(sm,k); it; ++it) {
				file_connection << it.row() << ", "; // row index
				file_connection << it.col() << ", "; // col index
				file_connection << it.value() << std::endl;
			}
		}
	} else {
		std::cout << "mat.innerSize: " << sm.innerSize() << "\n";
		std::cout << "mat.outerSize: " << sm.outerSize() << "\n";
		std::cout << "mat.nonZeros: " << sm.nonZeros() << "\n";
 	
	
		double* Values = sm.valuePtr();		  	// Pointer to the values
		int* InnerIndices = sm.innerIndexPtr();	// Pointer to the indices.
		int* OuterStarts = sm.outerIndexPtr();		// Pointer to the beginning of each inner vector
		// int* InnerNNZs = sm.innerNonZeroPtr();	// Not needed for compressed case
		
		
		
		for(int i = 0; i < sm.nonZeros(); ++i){
			file_connection << *(Values+i) << " ";
		}
		file_connection << std::endl;
		for(int i = 0; i < sm.nonZeros(); ++i){
			file_connection << *(InnerIndices+i) << " ";
		}
		file_connection << std::endl;
		for(int i = 0; i <= sm.outerSize(); ++i){
			file_connection << *(OuterStarts+i) << " ";
		}
		/*
		if(if_nnz){
			std::cout << std::endl;
			for(int i = 0; i < sm.innerSize(); ++i){
				std::cout << (*(InnerNNZs+i)) << " ";
			}
		}
		*/
	}
	
	
	
	file_connection << "\n";
	file_connection.close();
}









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
SpMat Hessian_mat(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                  const Vector_eig &TE, const Vector_eig &TR, 
                  const Vector_eig &sigma, const Matrix_eig_row &r, 
                  int n_x, int n_y, int n_z, MRF_param &MRF_obj, int with_MRF = 1, int verbose = 1){

	
	auto time_1_hess = std::chrono::high_resolution_clock::now();
	Debug1("Hessian calculation started");
	Matrix_eig_row v = v_mat(W, TE, TR);
	int n = n_x * n_y * n_z;
	int m = v.cols();
	double temp = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp31 = 0.0, tmp4 = 0.0;
	SpMat Gamma_inv;
	if(with_MRF){
		Gamma_inv = MRF_obj.Lambda(beta);
	}
	SpMat W_hess(3*n, 3*n);
	if(with_MRF){
		W_hess.reserve( VectorXi::Constant(3*n, 7*3) );
		// Reserve 7*3 non-zero's per column - https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
	} else {
		W_hess.reserve( VectorXi::Constant(3*n, 3) );
	}
	Debug1("Hessian matrix allocated");
	
	
	// First, the Kroneker prod term:
	if(with_MRF){
		SpMat Psi_inv_sp = Psi_inv.sparseView();
		W_hess = -Kron_Sparse_eig(Gamma_inv, Psi_inv_sp);
		//Debug1(" kron W_hess: \n");
		//show_head(MatrixXd(W_hess));
		// show_head_sp(W_hess);
		Debug0("MRF part done of Hessian!");
	}
	
	
	// Diagonal parts //
	int i = 0, i1 = 0, k = 0, k1 = 0, j = 0;
	Vector_eig temp_vec(3), temp_vec_1(3), temp_vec_2(3);;
	
	
	for(i = 0; i < n; ++i) {
	
		//if(i==100000 || i==300000 || i==500000 || i==700000 || i==900000 ){
		if(i % 10000 == 0){
			std::cout << "\n";
			Debug1("Hess matrix i: "<< i << ", j: " << j);
		}
		
		
		//temp_vec = W.row(i);
		for(k = 0; k < 3; ++k) {
			//for(k1 = 0; k1 < 3; ++k1) {
			for(k1 = k; k1 < 3; ++k1){
				
				temp = 0.;
				for(j = 0; j < m ; ++j) {
					
					tmp2 = r(i,j)/SQ(sigma(j));
					//if(i == 0){
						//Debug0(tmp2);
					//}
					tmp3 = - v(i,j)/SQ(sigma(j)) + tmp2 * besselI1_I0(tmp2 * v(i,j));
					//if(i == 0){
						//Debug0(tmp3);
					//}
					temp += tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1);
					//if(i == 0){
						//Debug0(simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1));  // problem
						//Debug0(tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1));
					//}
					
					
					//tmp2 *= v(i,j);
					//tmp3 = (1 + ratio_bessel_20(tmp2) - 2*SQ(besselI1_I0(tmp2)) );
					//tmp4 = -1/SQ(sigma(j)) +  0.5*SQ(r(i,j)/SQ(sigma(j)))*tmp3;
					// This is also valid
					
					
					
					
					// tmp4 = (-1)/SQ(sigma(j)) + SQ(tmp2) * h(tmp2*v(i, j)/SQ(sigma(j)));  // BUG
					tmp4 = (-1)/SQ(sigma(j)) + SQ(tmp2) * h(tmp2*v(i, j));
					// This is also valid
					
					//if(i == 0){
					//	Debug0(tmp2);
					//	Debug0(v(i,j));
					//	Debug0(SQ(sigma(j)));
					//	Debug0(tmp2*v(i, j)/SQ(sigma(j)));
					//	Debug0(h(tmp2*v(i, j)/SQ(sigma(j))));
					//	Debug0(tmp4);
					//}
					temp += tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
									simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1);
					//if(i == 0){
					//	Debug0(simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k));
					//	Debug0(simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
					//				simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1));
					//	Debug0(tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
					//				simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1));
					//	Debug0("Added:");
					//	Debug0(tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1));
					//	Debug0(tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
					//				simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1));
					//	Debug0((tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1) + 
					//				tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
					//					simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1)));
					//	Debug0("\n")
					//}
					//if(i == 0){
					//	Debug0(temp);
					//	Debug0("\n\n");
					//}
					
					if(k == k1 && temp > 0.1){
						Debug0("i: " << i << ", j: " << j << ", k: " << k);
						Debug0("W.row(i): " << W.row(i) << "\t r.row(i): " << r.row(i) << "\tsigma(j): " << sigma(j));
						Debug0("tmp3: " << tmp3 << "\t tmp4: " << tmp4);
						Debug0("Added: 1st part: " << tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1)
								<< ", 2nd part: " << tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
									simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1));
						Debug0("final: " << (tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1) + 
									tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
										simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1)) << "\n");
						if(temp > 100){
							Debug0("very high!!");
							// exit(EXIT_FAILURE);
							// change
						}
					}
					
				}
				// W_hess.insert(i+k*n, i+k1*n) = temp;		// old way - not very visually pleasing I guess.
				
				if(with_MRF){
					if(k == k1){
						W_hess.coeffRef(3 * i + k, 3 * i + k) += temp;
						// BUG? check negativity and positivity
						// minus added for Gamma * Psi						
					} else {
						W_hess.coeffRef(3 * i + k, 3 * i + k1) += temp;
						W_hess.coeffRef(3 * i + k1, 3 * i + k) += temp;
					}
				} else {
					if(k == k1){
						W_hess.insert(3 * i + k, 3 * i + k) = temp;
					} else {
						W_hess.insert(3 * i + k, 3 * i + k1) = temp;
						W_hess.insert(3 * i + k1, 3 * i + k) = temp;
					}
				}
			}
		}
	}
	
	
	
	W_hess.makeCompressed();
	//Debug1("W_hess: \n");
	//show_head(MatrixXd(-W_hess));
	// show_head_sp(-W_hess);
	
	auto time_2_hess = std::chrono::high_resolution_clock::now();
	auto duration_hess = std::chrono::duration_cast<std::chrono::seconds>(time_2_hess - time_1_hess);
	Debug1("Time taken total loop: " << duration_hess.count() << " seconds\n");
	Debug0("Hessian calculated with MRF");
	
	//show_head(W_hess);
	
	// return W_hess;	//3nx3n
	return (-W_hess);	//3nx3n
}
// Check sign please





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
* 
* 
* grad is changed.
* 
*/

void v_grad(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
             const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, const Matrix_eig_row &r, 
             int n_x, int n_y, int n_z, int i, int j, SpVec &grad){


	//SpMat grad(3*n_x*n_y*n_z, 1);
	//grad.reserve(VectorXi::Constant(1, 3));
	
	//for(int i1 = 3*i; i1 < 3 * i + 3; ++i1){
	//	grad.insert(i1, 1) = simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, i1 % 3);
	//}

	// SpVec grad(3*n_x*n_y*n_z);
	// Allocating this large vector might be cumbersome
	grad.setZero();
	
	
	for(int i1 = 3*i; i1 < 3 * i + 3; ++i1){
		// grad(i1) = simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, i1 % 3);
		grad.insert(i1) = simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, i1 % 3);
	}
	// return grad;
}






	// Export the result to a file:
	//	saveAsBitmap(x, n, argv[1]);







int choose(int n, int r){
	int tmp = n;
	for(int i = 1; i < r; ++i){
		tmp *= (n-i);
		tmp /= (i+1);
	}
	return tmp;
}


// https://stackoverflow.com/a/9430993
Matrix_eig combi(int n, int r){

	std::vector<bool> v(n);
	std::fill(v.begin(), v.begin() + r, true);
	int m = choose(n, r);
	Matrix_eig tmp = Matrix_eig::Zero(m, r);
	
	int k1 = 0, k2 = 0;
	
	do {
		k2 = 0;
		for (int i = 0; i < n; ++i) {
			if (v[i]) {
				// std::cout << (i + 1) << " ";
				tmp(k1, k2) = i; // + 1;
				k2++;
			}
		}
		k1++;
	} while (std::prev_permutation(v.begin(), v.end()));
	
	return tmp;
}






#endif	/* MAIN_HEADER */

