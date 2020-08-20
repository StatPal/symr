/**
* 
* To compile:

g++ scheme_new_OSL_EM_16_GEM.cpp -o test -I /usr/include/eigen3 -O3

g++ ~/MRI/Headers/TRY_EIGEN_5_NEW/scheme_new_OSL_EM_16_GEM.cpp -o scheme_new_OSL_EM_16_GEM -I ~/MRI/Headers -O3 -std=c++11

g++ scheme_new_OSL_EM_16_GEM.cpp -o scheme_new_OSL_EM_16_GEM -I ../eigen-3.3.7 -O3



./test ../Read_Data/new_phantom.nii Dummy_sd.txt 0
./scheme_new_OSL_EM_16_GEM ../Read_Data/new_phantom.nii Dummy_sd.txt 0

./scheme_new_EM_16_GEM ../Read_Data/new_phantom.nii Dummy_sd.txt 0


OSL EM algorithm
Psi and beta are updated at each iteration for this file. 


Changes:
Selectively, Matrix_eig to Matrix_eig_row
and, Eigen::MatrixXd to Matrix_eig_row
TMatrix to TMatrix_row
* 
*/




#include "scheme_new_numerical.hpp"
#include "Read_files_2.hpp"
#include "Init_value_6_numerical.hpp"

#include "../CppNumericalSolvers/include/cppoptlib/meta.h"
#include "../CppNumericalSolvers/include/cppoptlib/boundedproblem.h"
#include "../CppNumericalSolvers/include/cppoptlib/solver/lbfgsbsolver.h"

#include <ctime>
#include <iomanip>










/*
Penalised NEGATIVE log likelihood -- to be minimised
Matrix sizes: nx3, 3x3, 3(2)x1, mx1, mx1, mx1, nxm, ...
No change for cholesky inside the function
*/
double l_star(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta,
              const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, const Matrix_eig_row &r, 
              int n_x, int n_y, int n_z, MRF_param &MRF_obj){

	Matrix_eig_row v = v_mat(W, TE, TR);						// Can be passed
	int m = v.cols(), n = v.rows();
	double tmp2 = 0.0, tmp3 = 0.0, tmp1 = 0.0;

	//Rice part://
	int i = 0, j = 0;
	long double likeli_sum = 0.0;
	for(i = 0; i < n; ++i) {
		for(j = 0; j < m; ++j) {
			tmp2 = r(i,j)/SQ(sigma(j));
			tmp3 = (SQ(r(i,j))+SQ(v(i,j)))/SQ(sigma(j));
			tmp1 = logBesselI0(tmp2*v(i,j));
			likeli_sum += (log(tmp2) + tmp1 - 0.5*tmp3) ;
		}
	}
	
	//MRF part://
	likeli_sum += MRF_obj.MRF_log_likeli(W, Psi_inv, beta);
	
	// assert( ! std::isnan(-likeli_sum) );		// There are still NaN at W matrix after one iteration.
	return -likeli_sum;
}





/*
NEGATIVE Q_OSL fn, w.r.t. parameter per voxel -- to be minimised
Matrix sizes: nx3, 3x3, 3(2)x1, mx1, mx1, mx1, nxm, ...
*/
double Q_OSL_per_voxel(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                       const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, 
                       const Matrix_eig_row &r, const Matrix_eig_row &W_old, const Vector_eig &c_i,
                       int n_x, int n_y, int n_z, int i){

	
	Vector_eig v_i = Bloch_vec(W.row(i), TE, TR);
	Vector_eig v_old_i = Bloch_vec(W_old.row(i), TE, TR);
	int m = TE.size();
	double likeli_sum = 0.0, tmp2 = 0.0, tmp3 = 0.0;
	
	//Rice part://
	for(int j = 0; j < m; ++j) {
		tmp2 = r(i,j)*v_old_i(j)/SQ(sigma(j));
		tmp3 = besselI1_I0(tmp2);
		likeli_sum += v_i(j)*(- 0.5*v_i(j) + r(i,j)*tmp3)/SQ(sigma(j));
	}
	
	//MRF part://
	for(int k = 0; k < 3; ++k){
		likeli_sum -= c_i(k) * W(i, k);
	}
	
	//assert( ! std::isnan(-likeli_sum) );			// Don't assert, sobar mongol...
	return (-likeli_sum);
	
}




/*
* Negative Gradient of Penalised Q function per voxel - grad of J evaluated at old parameter
* Check simple_dee_v_ij_dee_W_ik again.
*/
Vector_eig Q_OSL_grad_per_voxel(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                                const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, 
                                const Matrix_eig_row &r, const Matrix_eig_row &W_old, const Vector_eig &c_i,
                                int n_x, int n_y, int n_z, int i){

	
	int m = TE.size();
	double temp = 0.0, tmp2 = 0.0, tmp3 = 0.0;
	Vector_eig W_grad(3);
	
	Vector_eig v_i = Bloch_vec(W.row(i), TE, TR);
	Vector_eig v_old_i = Bloch_vec(W_old.row(i), TE, TR);
		
	// Likelihood part
	for(int k = 0; k < 3; ++k){
		temp = 0.;
		for(int j = 0; j < m ; ++j){
			tmp2 = r(i,j)/SQ(sigma(j));
			tmp3 = -v_i(j)/SQ(sigma(j)) + tmp2 * besselI1_I0(tmp2*v_old_i(j));
			temp += tmp3 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k);
		}
		W_grad(k) = temp - c_i(k);
	}
	return (-W_grad);
	
}





/*
Penalised NEGATIVE Q fn, w.r.t. parameter of MRF -- to be minimised
* All parameters are not needed - keep them for uniformity?
*/
double Q_star_other_param(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta,
                          //const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, 
                          //const Matrix_eig_row &r, const Matrix_eig_row &W_old,
                          int n_x, int n_y, int n_z, MRF_param &MRF_obj){

	double likeli_sum = MRF_obj.MRF_log_likeli(W, Psi_inv, beta);
	return -likeli_sum;
}



/*
* Negative Gradient of Penalised Q function w.r.t. other parameters
*/
Vector_eig Q_grad_vec_other_parameter(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta,
                                      //const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, 
                                      //const Matrix_eig_row &r, const Matrix_eig_row &W_old,
                                      int n_x, int n_y, int n_z, MRF_param &MRF_obj){

	Vector_eig grad = MRF_obj.MRF_log_likeli_grad(W, Psi_inv, beta);	
	return (-grad);
}





/*
* Optim template for rows of W using partial fn:
*/
// Class definition with template
template<typename T>
class MRF_optim : public cppoptlib::BoundedProblem<T> {		// I guess it inherits
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;	 // Inherit the Vector typedef
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	typedef Matrix_eig_row TMatrix_row; 
	TMatrix_row r;
	TVector r2;
	MRF_param MRF_obj_optim;
	


// I guess, it's the declaration
  public:	
	// https://stackoverflow.com/a/18971392 :
	// Also see https://www.geeksforgeeks.org/when-do-we-use-initializer-list-in-c/
	MRF_optim(const TVector y_, const MRF_param &MRF_obj_optim) : 
		cppoptlib::BoundedProblem<T>(y_.size()),  
		r2(y_), 
		MRF_obj_optim(MRF_obj_optim) {}

	// I guess W can be included here if MRF is optimized for one time only.
	
	// OR, Let's do it with their style:
	// 
	// MRF_optim(const TVector y_, const MRF_param& MRF_obj_optim_) 
	// 	: cppoptlib::BoundedProblem<T>(y_.size()), r2(y_), MRF_obj_optim(MRF_obj_optim_){}


	int n_x, n_y, n_z;
	double beta_z = 0.1;												//Subrata - or get the value. 
	TVector lb, ub;
	TMatrix_row W;															// W here creating problem in optimization?
		
	
	
	
	
	// Track the best:
	Eigen::VectorXd current_best_param;
	double current_best_val = 1.0e+15;


	
	
	// x(8) has the following format:
	// First 6 elements for Cholesky decomposition of Psi_inv!
	// Then 2 element for beta_x and beta_y!
	T value(const TVector &x) {
	
		//Psi:
		Vector_eig temp_L = x.segment(0, 6);
		Matrix3d_eig L_mat = to_L_mat(temp_L);				// 6 element vec to Lower traingular mat: [a0, a1, a2; 0, a3, a4; 0, 0, a5]
		Matrix3d_eig Psi_inv = from_Cholesky(L_mat);		// Psi from Cholesky decomposition!
		
		//beta:
		Vector_eig beta = Vector_eig::Zero(3);
		beta(0) = x(6); beta(1) = x(7); beta(2) = beta_z;
		
		
		// check_bounds_vec(x, lb, ub);
		double fx = Q_star_other_param(W, Psi_inv, beta, n_x, n_y, n_z, MRF_obj_optim);
		//Debug2("- Q fn:" << fx);
		Debug2("x: " << std::setprecision(6) << x.transpose() << " \t& - Q fn:" << fx);
		
		// Track the best
		if(fx < current_best_val){
			current_best_param = x;
			current_best_val = fx;
		}
		return (fx);
	}

// Comment this Gradient part if you don't want to feed the gradient:

	
	
	// x(8) has the following format:
	// First 6 elements for Cholesky decomposition of Psi_inv!
	// Then 2 element for beta_x and beta_y!
	void gradient(const TVector &x, TVector &grad) {
		
		//Psi:
		Vector_eig temp_L = x.segment(0, 6);
		Matrix3d_eig L_mat = to_L_mat(temp_L);
		Matrix3d_eig Psi_inv = from_Cholesky(L_mat);
		
		//beta:
		Vector_eig beta = Vector_eig::Zero(3);
		beta(0) = x(6); beta(1) = x(7); beta(2) = beta_z;
		
		// check_bounds_vec(x, lb, ub);
		//Debug3("x in grad: "<< x.transpose());
		//Debug3("lb in grad: "<< lb.transpose());
		//Debug3("ub in grad: "<< ub.transpose());
		
		grad.noalias() = Q_grad_vec_other_parameter(W, Psi_inv, beta, n_x, n_y, n_z, MRF_obj_optim);
		
		
		
		// Chain rule for Cholesky:
		Vector_eig chain = grad.segment(0, 6);
		grad.segment(0, 6) = to_grad_Cholesky(temp_L)*chain;		//check transpose
		
		//Debug2("- grad:" << grad.transpose());
	}

};








/*
* Optim template for rows of W using partial fn:
*/
template<typename T>
class Likeli_optim : public cppoptlib::BoundedProblem<T> {			// Likeli_optim is inheriting fn from cppoptlib::BoundedProblem
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	typedef Matrix_eig_row TMatrix_row;
	TMatrix_row r;
	TVector r2;


  public:
	Likeli_optim(const TVector y_) : cppoptlib::BoundedProblem<T>(y_.size()), r2(y_){}





	int i, n_x, n_y, n_z;
	double beta_z = 0.1;												//Subrata - or get the value. 
	TVector TE, TR, sigma, beta, lb, ub, c_i;							// lb, ub are for extra check
	Matrix3d_eig Psi_inv;
	TMatrix_row W, W_old;													// W here creating problem in optimization?
	
	
	// Track the best:
	Eigen::VectorXd current_best_param;
	double current_best_val = 1.0e+15;


	T value(const TVector &x) {
		W.row(i) = x.transpose();
		//check_bounds_vec(x, lb, ub);
		double fx = Q_OSL_per_voxel(W, Psi_inv, beta, TE, TR, sigma, r, W_old, c_i, n_x, n_y, n_z, i);
		Debug2("x: " << x.transpose() << " \t& - Q fn:" << fx);
		
		// Track the best
		if(fx < current_best_val){
			current_best_param = x;
			current_best_val = fx;
		}
		return (fx);
	}

// Comment this Gradient part if you don't want to feed the gradient:

	void gradient(const TVector &x, TVector &grad) {
		W.row(i) = x;
		grad = Q_OSL_grad_per_voxel(W, Psi_inv, beta, TE, TR, sigma, r, W_old, c_i, n_x, n_y, n_z, i);
		Debug2("grad: " << grad.transpose() << "\n" );
	}

};












/*
* Main fn, currently one iteration is done. Change that with while loop
* Stopping criteria might seem confusing at first: 
* W_old is used to compare between new and previous iteration parameters
* and updated after each iteration
* whereas f.W_old is updated at each voxel update.
*/
void OSL_optim(Matrix_eig_row &W_init, Matrix3d_eig &Psi_inv, Vector_eig &beta, 
               const Vector_eig &TE_example, const Vector_eig &TR_example, 
               const Vector_eig &sigma, const Matrix_eig_row &r, 
               int n_x, int n_y, int n_z, double TE_scale, double TR_scale, 
               MRF_param &MRF_obj,
               int maxiter = 10, double abs_diff = 1e-6, int verbose = 0) {
// Change


	
	
	double old_val = 1.0e+15, old_likeli = 1.0e+15, current_best_likeli = 1.0e+15;
	int bad_count_o = 0, bad_count_o_2 = 0, bad_bound_1 = 0, bad_bound_2 = 0, nan_count = 0;
	
	
	///** First estimate other MRF parameters **///
	
	auto time_1_likeli = std::chrono::high_resolution_clock::now();
	MRF_optim<double> f_2(Eigen::VectorXd::Ones(8), MRF_obj);
	cppoptlib::LbfgsbSolver<MRF_optim<double>> solver_2;
	
	
	
	// *MRF based initial values:* //
	
	Eigen::VectorXd x_MRF(8), lb_MRF(8), ub_MRF(8), x_MRF_old(8);
	
	Matrix3d_eig L( Psi_inv.llt().matrixL() );
	x_MRF.segment(0, 6) = from_L_mat(L);
	x_MRF(6) = beta(0); x_MRF(7) = beta(1);
	x_MRF_old.noalias() = x_MRF;
	
	// Bound on Cholesky etc: (Why 255 btw? - Forgot - ad hoc I guess)
	ub_MRF.segment(0, 6+1+1) =  255*Vector_eig::Ones(6+1+1);
	lb_MRF.segment(0, 6+1+1) = -255*Vector_eig::Ones(6+1+1);
	lb_MRF(0) = 1e-5; lb_MRF(0+3) = 1e-5; lb_MRF(0+5) = 1e-5;	// Diagonals of Cholesky to be BDD away from 0
	lb_MRF(0+6) = 1e-5; lb_MRF(0+7) = 1e-5;						// beta forecefully made BDD away from 0
	
	//ub.segment(0, 6+1+1) *= 15;		// Get's sqaured
	//lb.segment(0, 6+1+1) = 0.00001*Vector_eig::Ones(6);
	
	f_2.setLowerBound(lb_MRF);	f_2.setUpperBound(ub_MRF);
	f_2.lb.noalias() = lb_MRF;	f_2.ub.noalias() = ub_MRF;		// Extra checks
	
	f_2.n_x = n_x; f_2.n_y = n_y; f_2.n_z = n_z;
	f_2.W.noalias() = W_init;
	
	
	
	
	// Subrata - Setting the parameters: new  -- (see simple_withoptions.cpp)
	cppoptlib::Criteria<double> crit_MRF = cppoptlib::Criteria<double>::defaults();
	crit_MRF.iterations = 1000;														// number of allowed iterations
	solver_2.setStopCriteria(crit_MRF);
	// Change 









	
	// * Voxel based initial values * //
	
	int iter = 0;
	Likeli_optim<double> f(Eigen::VectorXd::Ones(3));
	cppoptlib::LbfgsbSolver<Likeli_optim<double>> solver;			// For MRF parameters!
	
	// * Voxel based initial values * //
	
	Eigen::VectorXd x(3), lb(3), ub(3);
	
	//Bounds of rho, W1, W2:
	lb << 0.0, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	f.setLowerBound(lb);
	f.setUpperBound(ub);
	f.lb.noalias() = lb;	f.ub.noalias() = ub;								// Extra checks
	Debug2("lb: " << lb.transpose());
	Debug2("ub: " << ub.transpose());
	
	
	f.n_x = n_x; f.n_y = n_y; f.n_z = n_z;
	f.beta.noalias() = beta;
	f.Psi_inv.noalias() = Psi_inv;
	f.sigma.noalias() = sigma;	f.r.noalias() = r;	f.TE.noalias() = TE_example;	f.TR.noalias() = TR_example;
	f.W = W_init;
	Matrix_eig_row W_old = W_init;
	f.W_old.noalias() = W_old;
	
	
	
	// Subrata - Setting the parameters: new  -- (see simple_withoptions.cpp)
	cppoptlib::Criteria<double> crit_voxel = cppoptlib::Criteria<double>::defaults(); 	// Create a Criteria class to set the solver's stop conditions
	crit_voxel.iterations = 1000;														// Change the number of allowed iterations
	solver.setStopCriteria(crit_voxel);
	// Change: 
	
	
	
	
	
	
	
	// ** OSL-EM loop ** //
	
	while(iter < maxiter){
	
		
		Debug1("\n" << std::string(75, '-') << "\nIteration: " << iter++ << "\n");
		auto time_2_likeli = std::chrono::high_resolution_clock::now();
		
		
		
		// Track the best:
		f_2.current_best_val = 1.0e+15;
		
		//Print initial values:
		Debug2 ("x_MRF at first: " << x_MRF.transpose());
		Debug3 ("lb_MRF: " << lb_MRF.transpose());
		Debug3 ("ub_MRF: " << ub_MRF.transpose());
		Debug2 ("f(x) at first:");
		old_val = f_2.value(x_MRF_old);
		
		
		
		//Solve:
		solver_2.minimize(f_2, x_MRF);
		Debug2("argmin: " << x_MRF.transpose() << ";\tf(x) in argmin:");
		f_2(x_MRF);
		Debug2("Solver status: " << solver_2.status());
		Debug2("Final criteria values: " << "\n" << solver_2.criteria());
		Debug1("x_MRF: " << x_MRF.transpose());
		
	
		
		// Track the best: It's currently inside boundary
		x_MRF.noalias() = f_2.current_best_param;
		double fx_MRF = f_2.current_best_val;
		Debug2("best_param" << x_MRF.transpose() << "\t f(best_param): " << fx_MRF << 
				"\t old val:" << old_val << "\t diff: " << fx_MRF - old_val);
		
		if(fx_MRF >= old_val) {								//Compares best value inside
			Debug1("Value have not decreased!!\n" << " val: " << old_val << "; val: " << fx_MRF  << "\n");
			bad_count_o++;
			if(fx_MRF>old_val){
				bad_count_o_2++;
			}
		}
	
	
		// Calculated values: 
		// Could be shortened using f_2.value and f_2.gradient if gradient is not commented out: 
		// - but beta, Psi are needed for other optimization
		
		Vector_eig temp_L = x_MRF.segment(0, 6);
		Matrix3d_eig L_mat = to_L_mat(temp_L);
		Psi_inv.noalias() = from_Cholesky(L_mat);
		beta(0) = x_MRF(6); beta(1) = x_MRF(7); beta(2) = 0.1;
		
		/*
		double val_1 = Q_star_other_param(W_init, Psi_inv, beta, n_x, n_y, n_z);
		Debug2 ( "Optimized value: " << std::setprecision(15)  <<  val_1 << std::setprecision(6));
		
		Vector_eig grad1 = Q_grad_vec_other_parameter(W_init, Psi_inv, beta, n_x, n_y, n_z);
		// Chain rule for Cholesky:
		Vector_eig chain = grad1.segment(0, 6);
		grad1.segment(0, 6) = to_grad_Cholesky(temp_L)*chain;		//check transpose
		Debug2("- grad on optimized value:" << grad1.transpose());
		*/
		
		auto time_3_likeli = std::chrono::high_resolution_clock::now();
		auto duration_23 = std::chrono::duration_cast<std::chrono::seconds>(time_3_likeli - time_2_likeli);
		Debug1("Time taken for MRF optimization: " << duration_23.count() << " seconds\n");
		Debug1("MRF optimization ends!!");
		// * Optimization over other parameters ends * //
	






		// * Loop over voxels: * //
		
		
		// MRF contribution part:
		// One time per each loop - check removable or not:
		SpMat Gamma_inv = MRF_obj.Lambda(beta);
		Matrix_eig MRF_grad = Gamma_inv * W_old * Psi_inv;
		
		
		
		
		
		
		// Change: 
		for(int i = 0; i < r.rows()/1; ++i){
		//for(int i = 0; i < r.rows()/10000; ++i){
		//for(int i = 73; i < 75; ++i) {
		
			std::cout << "\n\n";
			if(i==100000 || i==200000 || i==300000 || i==400000 || i==500000 || i==600000 || i==700000 || i==800000 || i==900000 ){
				Debug1("i: "<< i);
			}
			
			f.i = i;
			f.c_i.noalias() = MRF_grad.row(i);
			x.noalias() = W_init.row(i);
			// check_bounds_vec(x, lb, ub);
			
			
			
			// Track the best:
			f.current_best_val = 1.0e+15;
			
			//Print initial values:
			Debug1 ("Value of i: " << i);
			Debug2 ("value of i: " << i << "\t x at first: " << x.transpose());
			Debug2 ("f(x) at first:");
			old_val = f.value(x);
			
			
			// Check derivative - new: 			// see Rosenbrock files
			/*
			bool probably_correct = f.checkGradient(x);
			if(probably_correct){
				Debug1(" Deriv is probably correct for voxel");
			} else {
				Debug1(" Deriv is probably NOT correct for voxel");
			}
			*/
			
			
			//Solve:
			solver.minimize(f, x);
			Debug1("argmin: " << x.transpose() << ";\tf(x) in argmin:");
			double fx = f(x);
			Debug2("Solver status: " << solver.status());	//Guess: bad reports: under constraints => grad is not ~0 
			Debug2("Final criteria values: " << "\n" << solver.criteria());
			
			
			// Track the best:
			// Don't do it now - best is going out of the box.(e.g., 1st elemet is going to ~6000)
			// x = f.current_best_param;
			// fx = f.current_best_val;
			Debug2("best_param: " << x.transpose() << "\t f(best_param): " << fx << 
					"\t old val: " << old_val << "\t diff: " << fx - old_val);
			
			
			if(fx >= old_val) {								//Compares best value inside
				Debug1("Value have not decreased!!\nold x:" << W_init.row(i) << " & val: " << old_val << 
						";\t x: " << x.transpose() << " val: " << fx << " i:" << i << "\n");
				bad_count_o++;
				if(fx > old_val){
					bad_count_o_2++;
				}
			} else {
				if(check_nan_vec(x) == 0){				// Added later, to catch NaN - Subrata
					W_init.row(i) = x;
				} else {
					nan_count++;
				}
			}
			
			
			// Restore values:
			f.W.row(i) = W_init.row(i);
			f.W_old.row(i) = W_init.row(i);
			
			
		}
		Debug0("Number of nan-voxels: " << nan_count << " at " << iter << "-th iter" );
		nan_count = 0;
		// * Voxel loop ends * //
		
		
		
		
		
		// *Checking stopping criterion with penalized negative log likelihood:* //
		int nan_count_2 = check_nan_W(W_init, W_old);					// Check this!
		Debug0("nan count in " << iter << "-th iteration: " << nan_count_2);
		current_best_likeli = l_star(W_init, Psi_inv, beta, TE_example, TR_example,
									 sigma, r, n_x, n_y, n_z, MRF_obj);
		
		
		if(current_best_likeli >= old_likeli){ 						// As everything is "-ve" log-likeli.
			Debug1("Value not decreased in EM loop!! old val: " << old_likeli << 
					";\t new val: " << current_best_likeli << " diff: " << current_best_likeli - old_likeli );
			//bad_count_o++;
		}
		Debug0(" Current likeli: " << -current_best_likeli << " Old likeli: " << -old_likeli << " diff: " << current_best_likeli - old_likeli );
		Debug1("Another iteration done\n\n");
		
		
		
		
		if(abs_sum(to_vector(W_old) - to_vector(W_init)) <= abs_diff && abs_sum(x_MRF_old - x_MRF) <= abs_diff){
		// remove the 2nd difference ???
			std::cout << "Stopped after " << iter << " iterations" << "\n";
			Debug1("W_old.row(73):" << W_old.row(73));
			Debug1("W_init.row(73):" << W_init.row(73));
			break;
		}
		Debug1("abs_sum(to_vector(W_old) - to_vector(W_init)):" << abs_sum(to_vector(W_old) - to_vector(W_init)));
		
		
		
		// Restore default values  ---- check other files also
		W_old.noalias() = W_init;
		x_MRF_old = x_MRF;
		old_likeli = current_best_likeli;
		
		
		auto time_4_likeli = std::chrono::high_resolution_clock::now();
		auto duration_34 = std::chrono::duration_cast<std::chrono::seconds>(time_4_likeli - time_3_likeli);
		Debug1("Time taken for one loop: " << duration_34.count() << " seconds\n");
	}
	if(iter > maxiter){
		Debug0("Max. iter reached for the ECM cycle");
	}
	// ** OSL-EM loop ends ** //
	
	
	
	
	
	
	std::cout << "\n\n";
	Debug0("Number of bad cases in Initial value determination:" << bad_count_o << 
			" and worse: " << bad_count_o_2 << 
			" and bad init bounds:" << bad_bound_1 << " and " << bad_bound_2);
	
	auto time_5_likeli = std::chrono::high_resolution_clock::now();
	auto duration_45 = std::chrono::duration_cast<std::chrono::seconds>(time_5_likeli - time_1_likeli);
	Debug1("Time taken for whole OSL-EM: " << duration_45.count() << " seconds\n");
	

}










int main(int argc, char * argv[]) {

	
	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);
	std::cout << "Current time: " << std::put_time(&tm, "%c %Z") << '\n';
	
	
	if (argc != 4) {
		fprintf(stderr, "\nUsage: %s <file_name> <SD_file_name> <will_write_to_a_file?> <temp_val> \n", argv[0]);
		exit(EXIT_FAILURE);
	}
	char *data_file, *sd_file;
	data_file = argv[1]; 	sd_file = argv[2]; 	char will_write = *(argv[3])-48;		// Converted from ascii
	short our_dim[8];
	
	
	


//	Vector_eig TE_example((Vector_eig(12) << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1).finished());
//	Vector_eig TR_example((Vector_eig(12) << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3).finished());




	Vector_eig TE_example((Vector_eig(18) << 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10).finished());
	Vector_eig TR_example((Vector_eig(18) << 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3).finished());
	double TE_scale = 1.01/TE_example.minCoeff();		// 1.01/0.03
	double TR_scale = 1.01/TR_example.minCoeff();		// 1.01/1.00
	Debug0("TE scale: " << TE_scale);
	Debug0("TR scale: " << TR_scale);
	TE_example *= TE_scale;
	TR_example *= TR_scale;
	//TE_scale, TR_scale are needed for determining the bounds
	
	
	
	Vector_eig lb(3), ub(3);
	lb << 0.0, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	Debug0("lb:" << lb.transpose());
	Debug0("ub:" << ub.transpose());
	
	
	
	double W1_init = exp(-1/(2.0*TR_scale));		// exp(-1/(2.0*1.01))
	double W2_init = exp(-1/(0.1*TE_scale));		// exp(-1/(0.1*1.01/0.03))
	





	
	Matrix_eig_row r = Preprocess_data(data_file, our_dim, will_write);
	Vector_eig sigma = read_sd(sd_file, our_dim[4]);
	Debug0("sigma: " << sigma.transpose());
	Debug2("Preprocessing done");
	

	
	// Divide into train and test:
	
	//std::vector<int> train_ind{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	//std::vector<int> test_ind{11};
	
	std::vector<int> train_ind{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	std::vector<int> test_ind{16, 17};
	
	
	Eigen::MatrixXd train(r.rows(), train_ind.size());
	Vector_eig TE_train(train_ind.size()), TR_train(train_ind.size()), sigma_train(train_ind.size());
	short our_dim_train[8];
	for(int i = 0; i < train_ind.size(); ++i) {
		train.col(i) = r.col(train_ind[i]);
		TE_train[i] = TE_example(train_ind[i]);
		TR_train[i] = TR_example(train_ind[i]);
		sigma_train[i] = sigma(train_ind[i]);
	}
	for(int i = 0; i < 8; ++i){
		our_dim_train[i] = our_dim[i];
	}
	our_dim_train[4] = (short)train_ind.size();		//our_dim[0] = 3 or 4
	
	Matrix_eig_row test(r.rows(), test_ind.size());
	Vector_eig TE_test(test_ind.size()), TR_test(test_ind.size()), sigma_test(test_ind.size());
	for(int i = 0; i < test_ind.size(); ++i){
		test.col(i) = r.col(test_ind[i]);
		TE_test[i] = TE_example(test_ind[i]);
		TR_test[i] = TR_example(test_ind[i]);
		sigma_test[i] = sigma(test_ind[i]);
	}
	
	
	
	// Temp results: Performance on the Init W: 
	Matrix_eig_row W_1st = Matrix_eig_row::Ones(our_dim_train[1]*our_dim_train[2]*our_dim_train[3], 3);
	W_1st.col(0) = train.rowwise().mean().transpose();
	// Knowing first non-zero index:
	int temp1 = 0;
	while(W_1st(temp1, 0) == 0.5){
		temp1++;
	}
	Debug1("First non-zero index: "<< temp1++);
	
	W_1st.col(1) *= W1_init;
	W_1st.col(2) *= W2_init;
	for(int i = 0; i < train.rows(); ++i){
		if(W_1st(i, 0) > 450){
			W_1st(i, 0) = 425.0;
		}
	}
	Matrix_eig perf_1 = Performance_test(W_1st, test, TE_test, TR_test, sigma_test, 1, 1);
	Matrix_eig perf_2 = Performance_test(W_1st, test, TE_test, TR_test, sigma_test, 3, 1);
	Matrix_eig perf_3 = Performance_test(W_1st, test, TE_test, TR_test, sigma_test, 1, 2);
	Matrix_eig perf_4 = Performance_test(W_1st, test, TE_test, TR_test, sigma_test, 3, 2);
	std::cout << "Performances over images: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images: " << perf_4.transpose() << "\n";
	
	
	
	
	
	// Least Sq: 
	// Change 
	int do_least_sq = 1;	// 0 Subrata -- least sq have better initial likelihood-but stucks and gives nan in some value
	Matrix_eig_row W_init = Init_val(train, TE_train, TR_train, our_dim_train, 
	                             TE_scale, TR_scale, W1_init, W2_init, do_least_sq, will_write);
	Debug1("W initial done");
	int nan_count_1st = check_nan_W(W_init, W_1st);
	Debug0("NAN count at first:" << nan_count_1st);
	show_head(W_init);
	perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
	std::cout << "Performances over images: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images: " << perf_4.transpose() << "\n";
	
	
	MRF_param MRF_obj_1(our_dim_train[1], our_dim_train[2], our_dim_train[3]);
	
	
	// Test:
	Matrix_eig_row W_LS = W_init;
	Debug1("abs diff between W's: " << abs_sum(to_vector(W_LS) - to_vector(W_init)));
	


	
	// Likelihood Based optimization:
	
	int n = our_dim[1]*our_dim[2]*our_dim[3];
	Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
	Vector_eig beta_init = 0.1*Vector_eig::Ones(3);						// beta_init	// Subrata multiply by 0.1
	OSL_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma, train, 
	          our_dim_train[1], our_dim_train[2], our_dim_train[3], TE_scale, TR_scale, MRF_obj_1);
	show_head(W_init);
	
	
	
	Debug1("abs diff between W's: " << abs_sum(to_vector(W_LS) - to_vector(W_init)));
		
	
	perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
	std::cout << "Performances over images: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images: " << perf_4.transpose() << "\n";
	
	
	
	std::time_t t2 = std::time(nullptr);
	std::tm tm2 = *std::localtime(&t2);
	std::cout << "Current time: " << std::put_time(&tm2, "%c %Z") << '\n';

	return 0;
}












