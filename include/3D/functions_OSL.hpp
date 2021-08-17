/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2021  Subrata Pal, Somak Dutta, Ranjan Maitra

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/


/* _OSL_HEADER_ */
#ifndef _OSL_HEADER_3D_
#define _OSL_HEADER_3D_




#include "../functions_gen.hpp"
#include "../read_files.hpp"
#include "../functions_LS_and_init_value.hpp"


#include "../../CppNumericalSolvers/include/cppoptlib/meta.h"
#include "../../CppNumericalSolvers/include/cppoptlib/boundedproblem.h"
#include "../../CppNumericalSolvers/include/cppoptlib/solver/lbfgsbsolver.h"

#include <ctime>
#include <iomanip>










/**
Penalised NEGATIVE log likelihood -- to be minimised
Matrix sizes: nx3, 3x3, 3(2)x1, mx1, mx1, mx1, nxm, ...
*/
double l_star_3D_OSL(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta,
              const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, const Matrix_eig_row &r, 
              MRF_param &MRF_obj, int penalized){

	Matrix_eig_row v = v_mat(W, TE, TR);						// Can be passed
	int m = v.cols(), n = v.rows();
	double tmp2 = 0.0, tmp3 = 0.0, tmp1 = 0.0;

	//Rice part://
	int i = 0, j = 0;
	long double likeli_sum = 0.0;
	for(i = 0; i < n; ++i) {
		for(j = 0; j < m; ++j) {
			tmp2 = r(i,j)/SQ(sigma(j));
			tmp3 = (SQ(r(i,j))+SQ(v(i,j)))/( 2*SQ(sigma(j)) );
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







/**
* Optim template for rows of W using partial fn:
*/
template<typename T>
class MRF_optim_3D_OSL : public cppoptlib::BoundedProblem<T> {		// I guess it inherits
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;	 // Inherit the Vector typedef
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	typedef Matrix_eig_row TMatrix_row;
	
	
	const TMatrix_row &W1;
	MRF_param &MRF_obj_optim;
	TMatrix tmp1, tmp2, tmp3;
	double fx;
	


  public:	
	MRF_optim_3D_OSL(const TMatrix_row &W1_, MRF_param &MRF_obj_optim_) : 
		cppoptlib::BoundedProblem<T>(1), 
		W1(W1_),
		MRF_obj_optim(MRF_obj_optim_), 
		tmp1(W1.transpose() * MRF_obj_optim_.H_1 * W1),
		tmp2(W1.transpose() * MRF_obj_optim_.H_2 * W1),
		tmp3(W1.transpose() * MRF_obj_optim_.H_3 * W1) {}
	

	TMatrix Psi_est, Psi_inv_est;
	TVector beta1 = (TVector(3) << 1, 1, 1).finished();
	
	
	void update_tmp(const TMatrix &W1){
		tmp1 = W1.transpose() * MRF_obj_optim.H_1 * W1;
		tmp2 = W1.transpose() * MRF_obj_optim.H_2 * W1;
		tmp3 = W1.transpose() * MRF_obj_optim.H_3 * W1;
	}
	

	// Get back the Psi_inv vector from beta vector
	TMatrix Psi_inv_mat(TVector &x) {
		beta1 << x(0), x(1), 1;
		Psi_est = (x(0) * tmp1 + x(1) * tmp2 + tmp3 )/MRF_obj_optim.n;
		Psi_inv_est = Psi_est.llt().solve(Matrix3d_eig::Identity(3, 3));
		return (Psi_inv_est);
	}
	
	// Objective function: 
	T value(const TVector &x) {
		beta1(0) = x(0); beta1(1) = x(1);
		Psi_est = (x(0) * tmp1 + x(1) * tmp2 + tmp3)/(MRF_obj_optim.n);
		fx = -(3 * MRF_obj_optim.sp_log_det_specific(beta1) - 
								MRF_obj_optim.n * log_det_3(Psi_est))/2;
		// Check the sign.
		return (fx);
	}
	
};








/**
* Optim template for rows of W using partial fn:
*/
template<typename T>
class Likeli_optim_3D_OSL : public cppoptlib::BoundedProblem<T> {
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	typedef Matrix_eig_row TMatrix_row;
	
	TMatrix_row r;
	
	TMatrix_row Theta;


  public:
	Likeli_optim_3D_OSL() : cppoptlib::BoundedProblem<T>(3){}


	int i;
	double beta_z = 1.0;
	TVector TE, TR, sigma, beta, lb, ub, c_i;								// lb, ub are for extra check
	Matrix3d_eig Psi_inv;
	TMatrix_row W, W_old;
	int penalized;
	
	
	// Track the best:
	Eigen::VectorXd current_best_param;
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
	


	
	Vector_eig v_i;
	
	void update_size(){
	 	v_i = Vector_eig::Zero(TE.size());
	}
	
	
	// Do not forget to update c_i for each i
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
			for(int k = 0; k < 3; ++k){
				likeli_sum -= c_i(k) * W(i, k);
			}
			// c_i = (Gamma_inv * W_old * Psi_inv).row(i)
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
		
		
		// Likelihood part: 
		for(int k = 0; k < 3; ++k){
			temp = 0.;
			for(int j = 0; j < m ; ++j){
				tmp2 = r(i,j)/SQ(sigma(j));
				tmp3 = -v_i(j)/SQ(sigma(j)) + tmp2 * Theta(i, j);
				temp += tmp3 * simple_dee_v_ij_dee_W_ik(x, TE, TR, j, k);
			}
			if(penalized){
				grad(k) = temp - c_i(k);
			} else {
				grad(k) = temp;
			}
		}
		
		grad = -grad;
		
		Debug2("grad: " << grad.transpose() << "\n" );
	}

};












/**
* The function for One Step Late estimation:  
* Inputs: 
	W_init:		W matrix, passed 
	Psi_inv:	Psi_inv matrix, passed
	beta: 		beta vector, passed
	TE_example: TE values for the train data
	TR_example: TR values for the train data
	sigma: 		sigma values for the train data
	r: 			Observed values for the pixels, n x m matrix
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
void OSL_optim_3D(Matrix_eig_row &W_init, Matrix3d_eig &Psi_inv, Vector_eig &beta, 
               const Vector_eig &TE_example, const Vector_eig &TR_example, 
               const Vector_eig &sigma, const Matrix_eig_row &r, 
               double r_scale, double TE_scale, double TR_scale, 
               MRF_param &MRF_obj,
               const Eigen::Matrix<char, Eigen::Dynamic, 1> &black_list, 
               int maxiter = 20, int penalized = 1, 
               double abs_diff = 0.1, double rel_diff = 1e-4, int verbose = 0, int verbose2 = 0) {
// Change

	if(verbose)
		std::cout << "\n\n\n";
	if(penalized){
		Debug0("Doing OSL-EM Estimate!");
	} else {
		Debug0("Doing EM Estimate!");
	}
	

	
	double old_val = 1.0e+15, old_likeli = 1.0e+15, current_best_likeli = 1.0e+15, fx = 0.0;
	int bad_count_o = 0, bad_count_o_2 = 0, bad_bound_1 = 0, bad_bound_2 = 0, nan_count = 0; 
	int n = r.rows(), m = r.cols();
	
//	Eigen::Matrix<char, Eigen::Dynamic, 1> black_list = Eigen::Matrix<char, Eigen::Dynamic, 1>::Ones(n);
//	
//	for(int i = 0; i < n; ++i){
//		for(int j = 0; j < m; ++j){
//			if(r(i, j) > 50){
//				black_list(i) = 0;
//				break;
//			}
//		}
//	}
//	Debug0("Number of possible background voxels: " << (black_list.sum()));
	
	
	
	
	///** First estimate other MRF parameters **///
	
	auto time_1_likeli = std::chrono::high_resolution_clock::now();
	//if(penalized){
		
		MRF_optim_3D_OSL<double> f_2(W_init, MRF_obj);
		cppoptlib::LbfgsbSolver<MRF_optim_3D_OSL<double>> solver_2;
		
		
		// *MRF based initial values:* //
		Vector_eig x_MRF(2), lb_MRF(2), ub_MRF(2), x_MRF_old(2);
		lb_MRF = Vector_eig::Constant(2, 1e-5); ub_MRF = Vector_eig::Constant(2, 1e+5);
		f_2.setLowerBound(lb_MRF);	f_2.setUpperBound(ub_MRF);
		x_MRF(0) = beta(0); x_MRF(1) = beta(1);
		x_MRF_old.noalias() = x_MRF;
		
		cppoptlib::Criteria<double> crit_MRF = cppoptlib::Criteria<double>::defaults();
		crit_MRF.iterations = 25;
		solver_2.setStopCriteria(crit_MRF);
		//Change 
		
	//}







	
	// * Voxel based initial values * //
	
	int iter = 0;
	Likeli_optim_3D_OSL<double> f;
	cppoptlib::LbfgsbSolver<Likeli_optim_3D_OSL<double>> solver;			// For MRF parameters!
	
	Eigen::VectorXd x(3), lb(3), ub(3);
	//Bounds of rho, W1, W2:
	lb << 0.0001, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	for(int i = 1; i < 3; ++i){
		if(lb[i]<1.0e-8){
			lb[i] = 1.0e-8;
		}
	}
	
	f.setLowerBound(lb);	f.setUpperBound(ub);
	f.lb.noalias() = lb;	f.ub.noalias() = ub;								// Extra checks
	Debug2("lb: " << lb.transpose());
	Debug2("ub: " << ub.transpose());
	
	
	f.update_penalized(penalized);
	f.beta.noalias() = beta;
	f.Psi_inv.noalias() = Psi_inv;
	f.sigma.noalias() = sigma;	f.r.noalias() = r;	f.TE.noalias() = TE_example;	f.TR.noalias() = TR_example;
	f.update_size();
	f.W.noalias() = W_init;
	Matrix_eig_row W_old = W_init;
	f.W_old.noalias() = W_old;
	
	f.Theta = Matrix_eig_row::Zero(n, m);
	// E_step: would give initial nonzero Theta
	f.E_step_update();
	
	
	
	
	// Subrata - Setting the parameters: new  -- (see simple_withoptions.cpp)
	cppoptlib::Criteria<double> crit_voxel = cppoptlib::Criteria<double>::defaults(); 	// Criteria class
	crit_voxel.iterations = 25;															// number of allowed iterations
	solver.setStopCriteria(crit_voxel);
	// Change maxiter 
	
	
	SpMat Gamma_inv = MRF_obj.Lambda(beta);
	Matrix_eig MRF_grad = Gamma_inv * W_old * Psi_inv;					// Subrata - change to Matrix_eig_row and check
	f.c_i = Vector_eig::Zero(3);		// Would be changed if penalized
	
	old_likeli = l_star_3D_OSL(W_init, Psi_inv, beta, TE_example, TR_example,
									 sigma, r, MRF_obj, penalized);
	
	
	
	// ** OSL-EM loop ** //
	
	while(iter < maxiter){
		
		iter++;
		if(verbose){
			Debug1("\n" << std::string(75, '-') << "\nIteration: " << iter << "\n");
		}
		
		
		
		if(penalized){
		
			f_2.update_tmp(W_init);
			
			
			//Print initial values:
			Debug2 ("x_MRF at first: " << x_MRF.transpose());
			Debug3 ("lb_MRF: " << lb_MRF.transpose());
			Debug3 ("ub_MRF: " << ub_MRF.transpose());
			Debug2 ("f(x) at first:");
			old_val = f_2.value(x_MRF_old);
			
		
		
		
			//Solve:
			solver_2.minimize(f_2, x_MRF);
			Debug2("argmin: " << x_MRF.transpose() << ";\tf(x) in argmin:");
			double fx_MRF = f_2(x_MRF);
			Debug2("Solver status: " << solver_2.status());
			Debug2("Final criteria values: " << "\n" << solver_2.criteria());
			if(verbose)
				Debug1("x_MRF: " << x_MRF.transpose());
			
			
			Debug2("best_param" << x_MRF.transpose() << "\t f(best_param): " << fx_MRF << 
					"\t old val:" << old_val << "\t diff: " << fx_MRF - old_val);
			if(fx_MRF >= old_val) {
				if(verbose){
					Debug1("Value have not decreased(MRF)!!\t" << " old val: " << old_val << "; new val: " << fx_MRF  << "\n");
				}
			bad_count_o++;
				if(fx_MRF>old_val){
					bad_count_o_2++;
				}
			}
			
			// Calculated values: 
			beta(0) = x_MRF(0); beta(1) = x_MRF(1); beta(2) = 1.0;
			Psi_inv = f_2.Psi_inv_mat(x_MRF);
			Debug0("MRF optimization done!");
			
			// * Optimization over other parameters ends * //
		
		}
		
		
		
		
		
		
		
		
		// * Loop over voxels: * //
		
		
		// MRF contribution part:
		// One time per each loop - check removable or not:
		if(penalized){
			Gamma_inv = MRF_obj.Lambda(beta);
			MRF_grad = Gamma_inv * W_old * Psi_inv;
		}
		
		
		// Change: 
		// Little bit different answer: - Subrata
		#pragma omp parallel for default(none) firstprivate(f, solver) private (x, old_val, fx) shared(W_init, bad_count_o, nan_count, bad_count_o_2, r, TE_example, TR_example, n, verbose, verbose2, black_list, penalized, MRF_grad, Rcpp::Rcout)
		for(int i = 0; i < n; ++i){
			if(i % 100000 == 0 ){
				if(verbose){
					//Rcpp::Rcout << std::endl;
					Debug1("i: "<< i);
				}
			}
			
			
			if(black_list(i) == 0){
				
				// Track the best:
				f.current_best_val = 1.0e+15;
				
				
				f.i = i;
				if(penalized) {
					f.c_i.noalias() = MRF_grad.row(i); 
				}
				x.noalias() = W_init.row(i);
				// check_bounds_vec(x, lb, ub);
				
				
				//Print initial values:
				Debug2 ("value of i: " << i << "\t x at first: " << x.transpose());
				Debug2 ("f(x) at first: ---------");
				old_val = f.value(x);
				
				
				// Check derivative - new: 			// see Rosenbrock files in the Optimization folder
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
				
				
				// Track the best:
				x = f.current_best_param;
				fx = f.current_best_val;
				
				
				Debug2("Solver status: " << solver.status());	//Guess: bad reports: under constraints => grad is not ~0 
				Debug2("Final criteria values: " << "\n" << solver.criteria());
				
				
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
				
				// Restore values:
				f.W.row(i) = W_init.row(i);
				f.W_old.row(i) = W_init.row(i);
				
				
			}
		}
		
		// E_step:
		f.E_step_update();
		std::cout << std::flush;
		
		
		
		
		if(nan_count>0){
			Debug0("Number of nan-voxels: " << nan_count << " at " << iter << "-th iter" );
		}
		nan_count = 0;
		if(verbose)
			Debug1("Voxel Loop ends!!");
		// * Voxel loop ends * //
		
		
		
		
		
		
		
		
		// *Checking stopping criterion* //
		
		// w.r.t. W 
		if(abs_sum(to_vector(W_old) - to_vector(W_init)) <= abs_diff){
			std::cout << "Stopped after " << iter << " iterations" << "\n";
			break;
		}
		if(verbose)
			Debug1("abs_sum(to_vector(W_old) - to_vector(W_init)):" << abs_sum(to_vector(W_old) - to_vector(W_init)));
		
		
		
		// with penalized negative log likelihood:
		current_best_likeli = l_star_3D_OSL(W_init, Psi_inv, beta, TE_example, TR_example,
									 sigma, r, MRF_obj, penalized);
		
		if(current_best_likeli >= old_likeli){ 						// As everything is "-ve" log-likeli.
			if(verbose){											// I guesss it is not good to have verbose here
				Debug1("Value not decreased in EM loop!! old val: " << old_likeli << 
						";\t new val: " << current_best_likeli << " diff: " << current_best_likeli - old_likeli);
			}
			//bad_count_o++;
		}
		if(verbose){
			Debug0(" Current likeli: " << -current_best_likeli << " Old likeli: " << -old_likeli 
							<< " diff: " << current_best_likeli - old_likeli );
			Debug0("rel. diff.: " << fabs(current_best_likeli - old_likeli)/fabs(current_best_likeli) << 
					"\t abs diff:" << fabs(current_best_likeli - old_likeli));
		}
		if(fabs(current_best_likeli - old_likeli)/fabs(current_best_likeli) <= rel_diff || iter == maxiter){
			std::cout << "Stopped after " << iter << " iterations (rel. diff.: " 
					<< fabs(current_best_likeli - old_likeli)/fabs(current_best_likeli) << ") abs diff:" 
					<< fabs(current_best_likeli - old_likeli) << "\n";
			break;
		}
		old_likeli = current_best_likeli;
		if(verbose)
			Debug1("Another iteration done\n\n");
		
		
		
		
		// Restore default values  ---- check other files also
		W_old.noalias() = W_init;
		
		
	}
	if(iter > maxiter){
		Debug0("Max. iter reached for the ECM cycle");
	}
	// ** OSL-EM loop ends ** //
	
	
	
	
	
	// ** Final Psi and beta ** //
	double tmp_sum = 0;
	tmp_sum = (beta(0) + beta(1) + beta(2)) * 2;
	beta /= tmp_sum;
	Psi_inv *= tmp_sum;
	
	
	// std::cout << "\n";
	if(verbose){
		Debug0("Number of bad cases in Initial value determination:" << bad_count_o << 
				" and worse: " << bad_count_o_2 << 
				" and bad init bounds:" << bad_bound_1 << " and " << bad_bound_2);
	}
	
	auto time_5_likeli = std::chrono::high_resolution_clock::now();
	auto duration_45 = std::chrono::duration_cast<std::chrono::seconds>(time_5_likeli - time_1_likeli);
	if(verbose){
		if(penalized){
			Debug1("Time taken for whole OSL-EM: " << duration_45.count() << " seconds\n");
		} else {
			Debug1("Time taken for whole EM: " << duration_45.count() << " seconds\n");
		}
	}
	
}





#endif /* !_OSL_HEADER_ */
