/* _AECM_HEADER_ */
#ifndef _AECM_HEADER_
#define _AECM_HEADER_




#include "functions_gen.hpp"
#include "read_files.hpp"
#include "functions_LS_and_init_value.hpp"


#include "../../CppNumericalSolvers/include/cppoptlib/meta.h"
#include "../../CppNumericalSolvers/include/cppoptlib/boundedproblem.h"
#include "../../CppNumericalSolvers/include/cppoptlib/solver/lbfgsbsolver.h"

#include <ctime>
#include <iomanip>








/*
Penalised NEGATIVE log likelihood -- to be minimised
*/
double l_star(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta,
              const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, const Matrix_eig_row &r, 
              int n_x, int n_y, int n_z, MRF_param &MRF_obj, int penalized){

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
	if(penalized){
		likeli_sum += MRF_obj.MRF_log_likeli(W, Psi_inv, beta);
	}
	
	return -likeli_sum;
}












/*
* Optim template for rows of W using partial fn:
*/
// Class definition with template
template<typename T>
class MRF_optim : public cppoptlib::BoundedProblem<T> {
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	typedef Matrix_eig_row TMatrix_row;
	
	
	const TMatrix_row &W1;
	MRF_param &MRF_obj_optim;
	TMatrix tmp1, tmp2;
	double fx;
	


  public:	
	MRF_optim(const TMatrix_row &W1_, MRF_param &MRF_obj_optim_) : 
		cppoptlib::BoundedProblem<T>(1), 
		W1(W1_),
		MRF_obj_optim(MRF_obj_optim_), 
		tmp1(W1.transpose() * MRF_obj_optim_.H_1 * W1),
		tmp2(W1.transpose() * MRF_obj_optim_.H_2 * W1) {}
	

	TMatrix Psi_est, Psi_inv_est;
	TVector beta1 = (TVector(3) << 1, 1, 0).finished();
	
	
	void update_tmp(const TMatrix &W1){
		tmp1 = W1.transpose() * MRF_obj_optim.H_1 * W1;
		tmp2 = W1.transpose() * MRF_obj_optim.H_2 * W1;
	}
	

	// Get back the Psi_inv vector from beta vector
	TMatrix Psi_inv_mat(TVector &x) {
		beta1 << x(0), 1, 0;
		Psi_est = (x(0) * tmp1 + tmp2 )/(MRF_obj_optim.n);
		Psi_inv_est = Psi_est.llt().solve(Matrix3d_eig::Identity(3, 3));
		return (Psi_inv_est);
	}
	
	// Objective function: 
	T value(const TVector &x) {
		beta1(0) = x(0);
		Psi_est = (x(0) * tmp1 + tmp2 )/(MRF_obj_optim.n);
		fx = -(3 * MRF_obj_optim.sp_log_det_specific(beta1) - 
								MRF_obj_optim.n * log_det_3(Psi_est))/2;
		// Check the sign.
		return (fx);
	}

};












/*
* Optim template for rows of W using partial fn:
*/
template<typename T>
class Likeli_optim : public cppoptlib::BoundedProblem<T> {
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	typedef Matrix_eig_row TMatrix_row;
	
	
	MRF_param &MRF_obj_optim;
	const TMatrix_row &r;
	
	TMatrix_row& Theta;
	
	TMatrix_row &W, &W_old;
	
	int i, n_x, n_y, n_z, penalized;
	double beta_z = 1.0;
	const TVector &lb, &ub, &sigma, &TE, &TR;						// lb, ub are for extra check
	TVector& beta;
	Matrix3d_eig& Psi_inv;
	TVector v_i;

	
	TVector MRF_grad = TVector::Zero(3);

	
  public:
	Likeli_optim(MRF_param &MRF_obj_optim_, const TMatrix_row& r_, TMatrix_row& Theta_, 
				 TMatrix_row& W_, TMatrix_row& W_old_,
				 int n_x_, int n_y_, int n_z_, int penalized_,
				 const TVector& lb_, const TVector& ub_, 
				 const TVector& sigma_, const TVector& TE_, const TVector& TR_,
				 TVector& beta_, Matrix3d_eig Psi_inv_) : 
		cppoptlib::BoundedProblem<T>(3), 
		MRF_obj_optim(MRF_obj_optim_), 
		r(r_), 
		Theta(Theta_), 
		W(W_), W_old(W_old_), 
		n_x(n_x_), n_y(n_y_), n_z(n_z_), penalized(penalized_),
		lb(lb_), ub(ub_),
		sigma(sigma_), TE(TE_), TR(TR_),
		beta(beta_),
		Psi_inv(Psi_inv_)
		{}
	// https://stackoverflow.com/questions/1828037/whats-the-point-of-g-wreorder
	// See this for reorder warning
	
	
	
	// Track the best:
	Vector_eig current_best_param;
	double current_best_val = 1.0e+15;
	
	
	void update_penalized(int val){
		penalized = val;
	}
	

	void E_step_update(){
		double tmp2 = 0.0;
		int m = TE.size(), n = W.rows();
		Vector_eig v_old_i = Vector_eig::Zero(m);
		for(int i = 0; i < n; ++i){
			Bloch_vec(W_old.row(i), TE, TR, v_old_i);
			for(int j = 0; j < m; ++j) {
				tmp2 = r(i,j)*v_old_i(j)/SQ(sigma(j));
				Theta(i, j) = besselI1_I0(tmp2);
			}
		}
	}

	
	
	
	void update_size(){
	 	v_i = Vector_eig::Zero(TE.size());
	}
	
	
	
	
	// MRF parts: 
	Matrix_eig tmp_i_1 = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_2 = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_Psi_inv_final = Matrix_eig::Zero(1, 3);
	Matrix_eig tmp_i_Psi_inv_new = Matrix_eig::Zero(1, 3);
	double tmp_i_coeff_1 = 0.0;
	
	void update_neighbours_likeli(int i){
	
		tmp_i_1 = Matrix_eig::Zero(1, 3);
		tmp_i_2 = Matrix_eig::Zero(1, 3);
		
		for (SpMat::InnerIterator it(MRF_obj_optim.H_1, i); it; ++it){
			if(it.row() != i){
				tmp_i_1 += W.row(it.row()) * it.value();		// * beta(0)
			}
		}
		tmp_i_1 = tmp_i_1 * beta(0);
		
		for (SpMat::InnerIterator it(MRF_obj_optim.H_2, i); it; ++it){
			if(it.row() != i){
				tmp_i_2 += W.row(it.row()) * it.value();
			}
		}
		tmp_i_2 = tmp_i_2;
		
		tmp_i_Psi_inv_final.noalias() = tmp_i_1 + tmp_i_2;
		tmp_i_Psi_inv_final = tmp_i_Psi_inv_final * Psi_inv;
		tmp_i_Psi_inv_final *= 2;
		
		tmp_i_coeff_1 = beta(0) * MRF_obj_optim.H_1.coeff(i, i) + 
						MRF_obj_optim.H_2.coeff(i, i);
	}
	
	double MRF_log_likeli_num_i_new(const Vector_eig &x) {
	
		tmp_i_Psi_inv_new.noalias() = (x.transpose() * Psi_inv);
		tmp_i_Psi_inv_new = tmp_i_Psi_inv_final + tmp_i_coeff_1 * tmp_i_Psi_inv_new;
		
		return ( -0.5 * (tmp_i_Psi_inv_new * x).value());
	}
	
	
	void MRF_grad_fn(const Vector_eig &x, Vector_eig &MRF_grad){
	
		MRF_grad = Psi_inv * x;
		MRF_grad *= tmp_i_coeff_1;
		MRF_grad += 0.5 * tmp_i_Psi_inv_final.transpose();
	}
	
	
	
	
	
	
	T value(const TVector &x) {
	
		W.row(i) = x.transpose();		
		Bloch_vec(W.row(i), TE, TR, v_i);
		int m = TE.size();
		double likeli_sum = 0.0;
		
		//Rice part://
		for(int j = 0; j < m; ++j) {
			likeli_sum += v_i(j)*(- 0.5*v_i(j) + r(i,j) * Theta(i, j))/SQ(sigma(j));
		}
		
		//MRF part://
		if(penalized){
			likeli_sum += MRF_log_likeli_num_i_new(x);
		}
		
		
		
		
		// Track the best:
		if((-likeli_sum) < current_best_val){
			if(check_bounds_vec_3(x, lb, ub) == 0){
				current_best_param = x;
				current_best_val = -likeli_sum;
			}
		}
		
		
		
		Debug2("x: " << x.transpose() << " \t&  Q fn:" << likeli_sum);
		return (-likeli_sum);
	}



// Comment this Gradient part if you don't want to feed the gradient:

	void gradient(const TVector &x, TVector &grad) {
	
		W.row(i) = x;
		
		
		int m = TE.size();
		double temp = 0.0, tmp2 = 0.0, tmp3 = 0.0;
		Bloch_vec(x, TE, TR, v_i);
		
	
		
		// MRF contribution part: 
		if(penalized){
			MRF_grad_fn(x, MRF_grad);
		}
		
		// Likelihood part: 
		for(int k = 0; k < 3; ++k){
			temp = 0.;
			for(int j = 0; j < m ; ++j){
				tmp2 = r(i,j)/SQ(sigma(j));
				tmp3 = -v_i(j)/SQ(sigma(j)) + tmp2 * Theta(i, j);
				temp += tmp3 * simple_dee_v_ij_dee_W_ik(x, TE, TR, j, k);
			}
			if(penalized){
				grad(k) = temp - MRF_grad(k);
			} else {
				grad(k) = temp;
			}
		}
		
		grad = -grad;
		Debug2("grad: " << grad.transpose() << "\n" );
	}
};












/*
* The function for One Step Late estimation:  
* Inputs: 
	W_init:		W matrix, passed 
	Psi_inv:	Psi_inv matrix, passed
	beta: 		beta vector, passed
	TE_example: TE values for the train data
	TR_example: TR values for the train data
	sigma: 		sigma values for the train data
	r: 			Observed values for the pixels, n x m matrix
	n_x, 
	n_y, 
	n_z: 		
	r_scale: 	scale for the r matrix, or equivalently rho.
	TE_scale: 	scale used for TE
	TR_scale: 	scale used for TR
	MRF_obj:	MRF_param object
	maxiter: 	Maximum number of iteration of EM algorithm - default value 20
	penalized: 	1(default) if penalization is used - 0 if not
	abs_diff: 	absolute difference bound for the EM algorithm
	rel_diff: 	relative difference bound for the EM algorithm for the log-likelihood
	verbose: 	verbose, default 0
	verbose2:	Secondary level of verbose - default 0
*/
void AECM_optim(Matrix_eig_row &W_init, Matrix3d_eig &Psi_inv, Vector_eig &beta, 
               const Vector_eig &TE_example, const Vector_eig &TR_example, 
               const Vector_eig &sigma, const Matrix_eig_row &r, 
               int n_x, int n_y, int n_z, double r_scale, double TE_scale, double TR_scale, 
               MRF_param &MRF_obj,
               int maxiter = 50, int penalized = 1, 
               double abs_diff = 1e-1, double rel_diff = 1e-5, int verbose = 0, int verbose2 = 0) {
// Change

	if(verbose)
		std::cout << "\n\n\n";
	if(penalized){
		Debug0("Doing AECM Estimate!");
	} else {
		Debug0("Doing EM Estimate!");
	}
	
	
	double old_val = 1.0e+15, old_likeli = 1.0e+15, current_best_likeli = 1.0e+15, fx = 0.0;
	int bad_count_o = 0, bad_count_o_2 = 0, bad_bound_1 = 0, bad_bound_2 = 0, nan_count = 0; 
	int n = r.rows(), m = r.cols();
	
	old_likeli = l_star(W_init, Psi_inv, beta, TE_example, TR_example,
							sigma, r, n_x, n_y, n_z, MRF_obj, penalized);
	
	
	Eigen::Matrix<char, Eigen::Dynamic, 1> black_list = Eigen::Matrix<char, Eigen::Dynamic, 1>::Ones(n);
	
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			if(r(i, j) > 50){
				black_list(i) = 0;
				break;
			}
		}
	}
	Debug0("Number of possible background voxels: " << (black_list.sum()));
	
	
	
	Eigen::Matrix<char, Eigen::Dynamic, 1> checkerboard_white = Eigen::Matrix<char, Eigen::Dynamic, 1>::Zero(n);
	int k = 0;
	for(int i = 0; i < MRF_obj.n_y_; ++i){
		for(int j = 0; j < MRF_obj.n_x_; ++j){					// Check the order
			checkerboard_white(k) = ((i % 2) + (j % 2)) % 2;
			k++;
		}
	}
	Debug0("Number of possible checkerboard white ones: " << (checkerboard_white.sum()));
	
	
	
	
	
	
	
	
	///** First estimate other MRF parameters **///
	
	auto time_1_likeli = std::chrono::high_resolution_clock::now();
	//if(penalized){
		
		MRF_optim<double> f_2(W_init, MRF_obj);
		cppoptlib::LbfgsbSolver<MRF_optim<double>> solver_2;
		
		
		// *MRF based initial values:* //
		Vector_eig x_MRF(1), lb_MRF(1), ub_MRF(1), x_MRF_old(1);
		lb_MRF(0) = 1e-5; ub_MRF(0) = 1e+5;	f_2.setLowerBound(lb_MRF);	f_2.setUpperBound(ub_MRF);
		x_MRF(0) = beta(0);
		x_MRF_old.noalias() = x_MRF;

		cppoptlib::Criteria<double> crit_MRF = cppoptlib::Criteria<double>::defaults();
		crit_MRF.iterations = 25;							// change
		solver_2.setStopCriteria(crit_MRF);		
	//}






	
	// * Voxel based initial values * //
	int iter = 0;
	Matrix_eig_row W_old = W_init;
	Matrix_eig_row W_old_reserve = W_old;
	Matrix_eig_row Theta = Matrix_eig_row::Zero(n, m);
	
	Vector_eig x(3), lb(3), ub(3);
	//Bounds of rho, W1, W2:
	lb << 0.0001, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	for(int i = 1; i < 3; ++i){
		if(lb[i] < 1.0e-8){
			lb[i] = 1.0e-8;
		}
	}
	
	
	
	
	
	
	
	
	// ** AECM loop ** //
	
	while(iter < maxiter){
		
		iter++;
		if(verbose){
			Debug1("\n" << std::string(75, '-') << "\nIteration: " << iter << "\n");
		}
		
		
		if(penalized){
		
			f_2.update_tmp(W_init);
			old_val = f_2.value(x_MRF_old);
			
			//Solve:
			solver_2.minimize(f_2, x_MRF);
			double fx_MRF = f_2.value(x_MRF);
			Debug2("Solver status: " << solver_2.status());
			Debug2("Final criteria values: " << "\n" << solver_2.criteria());
			if(verbose)
				Debug1("x_MRF: " << x_MRF.transpose());
			
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
			beta(0) = x_MRF(0); beta(1) = 1.0; beta(2) = 0.0;
			Psi_inv = f_2.Psi_inv_mat(x_MRF);
			Debug0("MRF optimization done!");
			// * Optimization over other parameters ends * //
		}
		
		
		
		
		// f.E_step_update();			// E_step: would give initial nonzero Theta
		
		double tmp2 = 0.0;
		Vector_eig v_old_i = Vector_eig::Zero(m);
		for(int i = 0; i < n; ++i){
			Bloch_vec(W_old.row(i), TE_example, TR_example, v_old_i);
			for(int j = 0; j < m; ++j) {
				tmp2 = r(i,j)*v_old_i(j)/SQ(sigma(j));
				Theta(i, j) = besselI1_I0(tmp2);
			}
		}
		
		
		
		
		// * Loop over voxels: * //
		#pragma omp parallel default(none) firstprivate(x, old_val, fx, v_old_i, tmp2) shared(MRF_obj, Theta, W_init, W_old, W_old_reserve,    r, n_x, n_y, n_z, penalized, lb, ub, sigma, TE_example, TR_example, Psi_inv, beta,   n, m, verbose, verbose2, bad_count_o, bad_count_o_2, nan_count, std::cout, checkerboard_white, black_list)		// Check v_old_i, tmp2 -- CAREFULLY - Subrata
		{
		
			// lb, ub, n_x, etc would be shared
			// beta, Psi_inv would be Private???? -- no, they are not changed -- shared
			Likeli_optim<double> f(MRF_obj, r, Theta, W_init, W_old, n_x, n_y, n_z, penalized, lb, ub,
									sigma, TE_example, TR_example, beta, Psi_inv);
			f.setLowerBound(lb);	f.setUpperBound(ub);
			f.update_size();
			
			
	
			cppoptlib::LbfgsbSolver<Likeli_optim<double>> solver;			// For MRF parameters!
			cppoptlib::Criteria<double> crit_voxel = cppoptlib::Criteria<double>::defaults();
			crit_voxel.iterations = 50;
			solver.setStopCriteria(crit_voxel);
			
			
			// * Checkerboard white loop: * //
			
			#pragma omp for
			for(int i = 0; i < n; ++i){
				if(i % 100000 == 0 ){
					if(verbose){
						std::cout << std::endl;
						Debug1("i: "<< i);
					}
				}
				
				if(checkerboard_white(i) == 1){
					if(black_list(i) == 0){
						
						// Track the best:
						f.current_best_val = 1.0e+15;
						
						
						f.i = i;
						if(penalized) {
							f.update_neighbours_likeli(i);
						}
						x.noalias() = W_init.row(i);
									
						//Print initial values:
						Debug2 ("value of i: " << i << "\t x at first: " << x.transpose());
						Debug2 ("f(x) at first:");
						old_val = f.value(x);
						
						// Check derivative:
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
						Debug2("argmin: " << x.transpose() << ";\tf(x) in argmin:");
						fx = f.value(x);
						Debug2("Solver status: " << solver.status());
						Debug2("Final criteria values: " << "\n" << solver.criteria());
						
						
						// Track the best:
						x = f.current_best_param;
						fx = f.current_best_val;
						
						Debug2("best_param: " << x.transpose() << "\t f(best_param): " << fx << 
								"\t old val: " << old_val << "\t diff: " << fx - old_val);
						
						
						if(fx >= old_val) {
							if(verbose2){
								Debug1("Value have not decreased!!\nold x:" << W_init.row(i) << " & val: " << old_val << 
										";\t x: " << x.transpose() << " val: " << fx << " i:" << i << "\n");
							}
							bad_count_o++;
							if(fx > old_val){
								bad_count_o_2++;
							}
						} else {
							if(check_nan_vec(x) == 0){
								W_init.row(i) = x;
							} else {
								Debug1("nan in EM estimate. \n" << "i: " << i << ", x: " << x.transpose() << 
										"\nr.row(i): " << r.row(i));
								nan_count++;
							}
						}
						
						// Restore values:	// I guess this is not necessary now - Check - Subrata
						// f.W.row(i) = W_init.row(i);
						f.W_old.row(i) = W_init.row(i);
					}
				}
			}
			
			
			// Qn: Update W_old?????? - Subrata - check 
			// W_old = W_init;
			
			
			
			#pragma omp barrier
			
			// E_step:
			#pragma omp for
			for(int i = 0; i < n; ++i){
				Bloch_vec(W_old.row(i), TE_example, TR_example, v_old_i);
				for(int j = 0; j < m; ++j) {
					tmp2 = r(i,j)*v_old_i(j)/SQ(sigma(j));
					Theta(i, j) = besselI1_I0(tmp2);
				}
			}
			
			std::cout << std::flush;
			// * Checkerboard white ends * //
			
			
			
			
			
			
			#pragma omp barrier
			
			// * Checkerboard black loop: * //
			#pragma omp for
			for(int i = 0; i < n; ++i){
				if(i % 100000 == 0 ){
					if(verbose){
						std::cout << std::endl;
						Debug1("i: "<< i);
					}
				}
				
				if(checkerboard_white(i) == 0){
					if(black_list(i) == 0){
						// Track the best:
						f.current_best_val = 1.0e+15;
						
						f.i = i;
						if(penalized) {
							f.update_neighbours_likeli(i);
						}
						x.noalias() = W_init.row(i);
						
									
						//Print initial values:
						Debug2 ("value of i: " << i << "\t x at first: " << x.transpose());
						Debug2 ("f(x) at first:");
						old_val = f.value(x);
						
						
						// Check derivative
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
						Debug2("argmin: " << x.transpose() << ";\tf(x) in argmin:");
						double fx = f.value(x);
						Debug2("Solver status: " << solver.status());
						Debug2("Final criteria values: " << "\n" << solver.criteria());
						
						
						// Track the best:
						x = f.current_best_param;
						fx = f.current_best_val;
						
						
						Debug2("best_param: " << x.transpose() << "\t f(best_param): " << fx << 
								"\t old val: " << old_val << "\t diff: " << fx - old_val);
						
						
						if(fx >= old_val) {								//Compares best value inside
							if(verbose2){
								Debug1("Value have not decreased!!\nold x:" << W_init.row(i) << " & val: " << old_val << 
										";\t x: " << x.transpose() << " val: " << fx << " i:" << i << "\n");					
							}
							bad_count_o++;
							if(fx > old_val){
								bad_count_o_2++;
							}
						} else {
							if(check_nan_vec(x) == 0){				// Added later, to catch NaN - Subrata
								W_init.row(i) = x;
							} else {
								Debug1("nan in EM estimate. \n" << "i: " << i << ", x: " << x.transpose() << 
										"\nr.row(i): " << r.row(i));
								nan_count++;
							}
						}
						
						// Restore values:	// I guess this is not necessary now - Check - Subrata
						// f.W.row(i) = W_init.row(i);
						f.W_old.row(i) = W_init.row(i);
					}
				}
			}
			
			
			
			#pragma omp barrier
			// E_step:
			#pragma omp for
			for(int i = 0; i < n; ++i){
				Bloch_vec(W_old.row(i), TE_example, TR_example, v_old_i);
				for(int j = 0; j < m; ++j) {
					tmp2 = r(i,j)*v_old_i(j)/SQ(sigma(j));
					Theta(i, j) = besselI1_I0(tmp2);
				}
			}
			
			std::cout << std::flush;
			// * Checkerboard black ends * //
			
		}
		
		
		
		
		
		
		
		
		if(nan_count > 0){
			Debug0("Number of nan-voxels: " << nan_count << " at " << iter << "-th iter" );
		}
		nan_count = 0;
		if(verbose)
			Debug1("Voxel Loop ends!!");
		// * Voxel loop ends * //
		
		
		
		
		
		// *Checking stopping criterion:* //
		
		// w.r.t. W 
		if(abs_sum(to_vector(W_old_reserve) - to_vector(W_init)) <= abs_diff){
			std::cout << "Stopped after " << iter << " iterations" << "\n";			// This is weird to stop at
			if(verbose)
				Debug1("abs diff. (W_old - W_new):" << abs_sum(to_vector(W_old_reserve) - to_vector(W_init)));
			break;
		}
		if(verbose)
			Debug1("abs diff. (W_old - W_new):" << abs_sum(to_vector(W_old_reserve) - to_vector(W_init)));
		
		
				
		// with penalized negative log likelihood:
		current_best_likeli = l_star(W_init, Psi_inv, beta, TE_example, TR_example,
									 sigma, r, n_x, n_y, n_z, MRF_obj, penalized);
		
		
		if(current_best_likeli >= old_likeli){ 						// As everything is "-ve" log-likeli.
			if(verbose){											// I guesss it is not good to have verbose here
				Debug1("Value not decreased in EM loop!! old val: " << old_likeli << 
						";\t new val: " << current_best_likeli << 
						" diff: " << current_best_likeli - old_likeli);				
			}
			//bad_count_o++;
		}
		if(verbose){
			Debug0("Current log likeli: " << -current_best_likeli << ", Old log likeli: " << -old_likeli 
					<< ", diff.: " << current_best_likeli - old_likeli  << 
					",  rel. diff.: " << fabs(current_best_likeli - old_likeli)/fabs(current_best_likeli));
		}
		if(fabs(current_best_likeli - old_likeli)/fabs(current_best_likeli) <= rel_diff || iter == maxiter){
			std::cout << "\nStopped after " << iter << " iterations (rel. diff.: " 
					<< fabs(current_best_likeli - old_likeli)/fabs(current_best_likeli) << ", abs diff: " 
					<< fabs(current_best_likeli - old_likeli) << ")\n";
			break;
		}
		old_likeli = current_best_likeli;
		if(verbose)
			Debug1("Another iteration done\n\n");
		
		
		
		
		
		// Restore default values  ---- check other files also
		W_old_reserve.noalias() = W_init;
		
		
	}
	if(iter > maxiter){
		Debug0("Max. iter reached for the ECM cycle");
	}
	// ** AECM loop ends ** //
	
	
	
	
	
	
	
	// std::cout << "\n";
	if(verbose){
		Debug1("Number of bad cases in optimization:" << bad_count_o << 
				" and worse: " << bad_count_o_2 << 
				" and bad init bounds:" << bad_bound_1 << " and " << bad_bound_2);		
	}
	
	auto time_5_likeli = std::chrono::high_resolution_clock::now();
	auto duration_45 = std::chrono::duration_cast<std::chrono::seconds>(time_5_likeli - time_1_likeli);
	if(verbose){
		if(penalized){
			Debug1("Time taken for whole AECM: " << duration_45.count() << " seconds\n");
		} else {
			Debug1("Time taken for whole EM: " << duration_45.count() << " seconds\n");
		}
	}
	
}




#endif /* !_AECM_HEADER_ */
