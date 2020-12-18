/**
* 
* OSL EM algorithm
* Psi and beta are updated at every loop
* E step is done seperately in that class of optimization. 


* 


* To compile:

g++ All_possible_26_GEM.cpp -o all_26_2D -I /usr/include/eigen3 -O3 -lgsl -lgslcblas -lm -fopenmp




./all_26_2D ../Read_Data/new_phantom.nii Dummy_sd.txt 0

./all_26_2D ../data/new_phantom.nii Dummy_sd.txt 0

./all_26_2D ../Read_Data/small_phantom.nii Dummy_sd.txt 0

nohup ./all_26_2D ../Read_Data/new_phantom.nii Dummy_sd.txt 0 > all_26_2D.out & 




Changes:

Black listed pixels


MRF estimation in a new way



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
#include <set>









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
	
	return -likeli_sum;
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
		Psi_est = (x(0) * tmp1 + tmp2 )/(MRF_obj_optim.n);				//  not n * 3
		Psi_inv_est = Psi_est.llt().solve(Matrix3d_eig::Identity(3, 3));
		return (Psi_inv_est);
	}
	
	// Objective function: 
	T value(const TVector &x) {
		beta1(0) = x(0);
		Psi_est = (x(0) * tmp1 + tmp2 )/(MRF_obj_optim.n);	// I guess there would be an additional 3. Check!
//		Psi_inv_est = Psi_est.llt().solve(Matrix3d_eig::Identity(3, 3));
//		double fx = -(3 * MRF_obj_optim.sp_log_det_specific(beta1) + 
//								MRF_obj_optim.n * log_det_3(Psi_inv_est))/2;
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
class Likeli_optim : public cppoptlib::BoundedProblem<T> {			// Likeli_optim is inheriting fn from cppoptlib::BoundedProblem
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	typedef Matrix_eig_row TMatrix_row;
	
	TMatrix_row r;
	
	TMatrix_row Theta;		// new


  public:
	Likeli_optim() : cppoptlib::BoundedProblem<T>(3){}


	int i, n_x, n_y, n_z;
	double beta_z = 0.0, beta_y = 1.0;										//Subrata - or get the value. 
	TVector TE, TR, sigma, beta, lb, ub, c_i;								// lb, ub are for extra check
	Matrix3d_eig Psi_inv;
	TMatrix_row W, W_old;													// W here creating problem in optimization?
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
		double likeli_sum = 0.0, tmp2 = 0.0, tmp3 = 0.0;
		
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












/*
* Main fn, 
* Stopping criteria might seem confusing at first: 
* W_old is used to compare between new and previous iteration parameters
* and updated after each iteration
* whereas f.W_old is updated at each voxel update.
*/
void OSL_optim(Matrix_eig_row &W_init, Matrix3d_eig &Psi_inv, Vector_eig &beta, 
               const Vector_eig &TE_example, const Vector_eig &TR_example, 
               const Vector_eig &sigma, const Matrix_eig_row &r, 
               int n_x, int n_y, int n_z, double r_scale, double TE_scale, double TR_scale, 
               MRF_param &MRF_obj,
               int maxiter = 20, int penalized = 1, 
               double abs_diff = 0.1, double rel_diff = 1e-4, int verbose = 0, int verbose2 = 0) {
// Change abs_diff, maxiter etc

	if(verbose)
		std::cout << "\n\n\n";
	if(penalized){
		Debug0("Doing OSL-EM Estimate!");
	} else {
		Debug0("Doing EM Estimate!");
	}
	
	
	double old_val = 1.0e+15, old_likeli = 1.0e+15, current_best_likeli = 1.0e+15, fx;
	int bad_count_o = 0, bad_count_o_2 = 0, bad_bound_1 = 0, bad_bound_2 = 0, nan_count = 0; 
	int n = r.rows(), m = r.cols();
	
	Eigen::Matrix<char, Dynamic, 1> black_list = Eigen::Matrix<char, Dynamic, 1>::Ones(n);
	
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			if(r(i, j) > 50){
				black_list(i) = 0;
				break;
			}
		}
	}
	Debug0("Number of possible background voxels: " << black_list.sum());
	
	
	
	
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
		crit_MRF.iterations = 50;
		solver_2.setStopCriteria(crit_MRF);
		
	//}







	
	// * Voxel based initial values * //
	
	int iter = 0;
	Likeli_optim<double> f;
	cppoptlib::LbfgsbSolver<Likeli_optim<double>> solver;			// For MRF parameters!
	
	// * Voxel based initial values * //
	
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
	
	
	f.n_x = n_x; f.n_y = n_y; f.n_z = n_z;
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
	crit_voxel.iterations = 80;															// number of allowed iterations
	solver.setStopCriteria(crit_voxel);
	// Change maxiter 
	
	
	SpMat Gamma_inv = MRF_obj.Lambda(beta);
	Matrix_eig MRF_grad = Gamma_inv * W_old * Psi_inv;					// Subrata - change to Matrix_eig_row and check
	f.c_i = Vector_eig::Zero(3);		// Would be changed if penalized
	
	old_likeli = l_star(W_init, Psi_inv, beta, TE_example, TR_example,
									 sigma, r, n_x, n_y, n_z, MRF_obj);
	
	
	
	// ** OSL-EM loop ** //
	
	while(iter < maxiter){
		
		iter++;
		if(verbose){
			Debug1("\n" << std::string(75, '-') << "\nIteration: " << iter << "\n");
		}
		auto time_2_likeli = std::chrono::high_resolution_clock::now();
		
		
		
		if(penalized){
		
			// f_2.W.noalias() = W_init;
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
			beta(0) = x_MRF(0); beta(1) = 1.0; beta(2) = 0.0;
			Psi_inv = f_2.Psi_inv_mat(x_MRF);
			Debug0("MRF optimization done!");
			
			// auto time_2_likeli = std::chrono::high_resolution_clock::now();
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
		#pragma omp parallel for default(none) firstprivate(f, solver, TE_example, TR_example) private (x, old_val, fx) shared(W_init, bad_count_o, nan_count, bad_count_o_2, r, n, verbose, verbose2, black_list, penalized, MRF_grad, std::cout)
		for(int i = 0; i < n; ++i){
		//for(int i = 73; i < 75; ++i) {
		
			if(i % 10000 == 0 ){
				if(verbose){
					std::cout << std::endl;
					Debug1("i: "<< i);
				}
			}
			
			
			if(black_list(i) == 0){
			
			
				auto time_1_voxel = std::chrono::high_resolution_clock::now();
				
				// Track the best:
				// double current_best_val = 1.0e+15;
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
				fx = f.value(x);
				
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
				
				
				auto time_2_voxel = std::chrono::high_resolution_clock::now();
				auto duration_12_voxel = std::chrono::duration_cast<std::chrono::microseconds>(time_2_voxel - time_1_voxel);
				//if(penalized)
				//	Debug1("Time taken for MRF part optim: " << duration_12_voxel.count() << " microseconds\n");
			}
			
		}
		
		// E_step:
		f.E_step_update();
		std::cout << std::flush;
		
		
		
		
		if(nan_count>0){
			Debug0("Number of nan-voxels: " << nan_count << " at " << iter << "-th iter" );
		}
		nan_count = 0;
		auto time_3_likeli = std::chrono::high_resolution_clock::now();
		//auto duration_23 = std::chrono::duration_cast<std::chrono::microseconds>(time_3_likeli - time_2_likeli);
		//Debug1("Time taken for 1 OSL-EM loop with " << r.rows() << " rows: " << duration_23.count() << " microseconds");
		if(verbose)
			Debug1("Voxel Loop ends!!");
		// * Voxel loop ends * //
		
		
		
		
		
		
		
		
		// *Checking stopping criterion* //
		
		// w.r.t. W 
		if(abs_sum(to_vector(W_old) - to_vector(W_init)) <= abs_diff){
			std::cout << "Stopped after " << iter << " iterations" << "\n";
			// Debug1("W_old.row(73):" << W_old.row(73));
			// Debug1("W_init.row(73):" << W_init.row(73));
			break;
		}
		if(verbose)
			Debug1("abs_sum(to_vector(W_old) - to_vector(W_init)):" << abs_sum(to_vector(W_old) - to_vector(W_init)));
		
		
		
		// with penalized negative log likelihood:
		current_best_likeli = l_star(W_init, Psi_inv, beta, TE_example, TR_example,
									 sigma, r, n_x, n_y, n_z, MRF_obj);
		
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
		
		
		
		
		
		
		
		auto time_4_likeli = std::chrono::high_resolution_clock::now();
		auto duration_34 = std::chrono::duration_cast<std::chrono::seconds>(time_4_likeli - time_3_likeli);
		if(verbose)
			Debug1("Time taken for MRF part optim: " << duration_34.count() << " seconds\n");
	}
	if(iter > maxiter){
		Debug0("Max. iter reached for the ECM cycle");
	}
	// ** OSL-EM loop ends ** //
	
	
	
	
	
	
	// std::cout << "\n";
	if(verbose){
		Debug0("Number of bad cases in Initial value determination:" << bad_count_o << 
				" and worse: " << bad_count_o_2 << 
				" and bad init bounds:" << bad_bound_1 << " and " << bad_bound_2);
	}
	
	auto time_5_likeli = std::chrono::high_resolution_clock::now();
	auto duration_45 = std::chrono::duration_cast<std::chrono::microseconds>(time_5_likeli - time_1_likeli);
	if(verbose){
                if(penalized){
                        Debug1("Time taken for whole OSL-EM: " << duration_45.count() << " microseconds\n");
                } else {
                        Debug1("Time taken for whole EM: " << duration_45.count() << " microseconds\n");
                }
        }	
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
	
	
	
	
	// Reading the data: 
	Matrix_eig_row r = Preprocess_data(data_file, our_dim, will_write);
	Vector_eig sigma = read_sd(sd_file, our_dim[4]);
	
	// Scaled: r, sigma, ub would change.
	double r_scale = r.maxCoeff();
	r_scale = 10.0;
	r.array() /= r_scale;
	sigma.array() /= r_scale;
	
	Debug0("sigma: " << sigma.transpose());
	Debug2("Preprocessing done");
	
	
	
	
	
	
	



	Vector_eig TE_example((Vector_eig(18) << 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10).finished());
	Vector_eig TR_example((Vector_eig(18) << 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3).finished());
	// 1.01 -> 2.01
	double TE_scale = 2.01/TE_example.minCoeff();		// 1.01/0.03
	double TR_scale = 2.01/TR_example.minCoeff();		// 1.01/1.00
	Debug0("r_scale: " << r_scale);
	Debug0("TE scale: " << TE_scale);
	Debug0("TR scale: " << TR_scale);
	TE_example *= TE_scale;
	TR_example *= TR_scale;
	//TE_scale, TR_scale are needed for determining the bounds
	
	
	
	Vector_eig lb(3), ub(3);
	lb << 0.0001, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	for(int i = 1; i < 3; ++i){
		if(lb[i]<1.0e-8){
			lb[i] = 1.0e-8;
		}
	}
	Debug0("lb:" << lb.transpose());
	Debug0("ub:" << ub.transpose());
	
	
	
	double W1_init = exp(-1/(2.0*TR_scale));		// exp(-1/(2.0*1.01))
	double W2_init = exp(-1/(0.1*TE_scale));		// exp(-1/(0.1*1.01/0.03))
	



	// To be stored in the file: 
	std::ofstream file_performance;
	file_performance.open ("result/Performances_26_all_3.csv");
	
	
	Matrix_eig perf_1, perf_2, perf_3, perf_4;

	// file_performance << "Test\n";
	std::cout << std::flush;
	
	
	
	int m_total = 18;
	std::vector<int> whole_ind = {};
	for(int i = 0; i < m_total; ++i)
		whole_ind.push_back(i);
	
	
	int r1 = m_total - 3;
//	int r1 = 1;
//	for(int r1 = 1; r1 < m_total - 3; ++r1){			// r1 is the test set size
		
		int m_choose_r1 = choose(m_total, r1);
		Debug0("m_choose_r: " << m_choose_r1);
		Matrix_eig tmp_combi = combi(m_total, r1);
		
		std::vector<int> train_ind = {};
		
		
		for(int i = 0; i < m_choose_r1; ++i){
		
			Debug1("i = " <<  i);
		
			// Divide into train and test:
			train_ind = {};
			
			// Test: 
			std::vector<int> test_ind = {};
			for(int j = 0; j < r1; ++j)
				test_ind.push_back( tmp_combi(i, j) );
			
			std::set_difference(whole_ind.begin(), whole_ind.end(), 
								test_ind.begin(), test_ind.end(),
                        		std::inserter(train_ind, train_ind.begin()));
			
			
			
			Matrix_eig_row train(r.rows(), train_ind.size());
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
			
			
			std::set<double> TE_unique_train{TE_train.data(), TE_train.data() + TE_train.size()};
			std::set<double> TR_unique_train{TR_train.data(), TR_train.data() + TR_train.size()};
			
			
//			if( ( TE_unique_train.size() >= 3 ) && (TR_unique_train.size() >= 3)){
			if( ( TE_unique_train.size() >= 2 ) && (TR_unique_train.size() >= 2)){
				
				std::cout << "\ntrain: \t";
				for(int n1 : train_ind) {
	    	 	   std::cout << n1 << ' ';
	    		}
				std::cout << "\ntest: \t";
				for(int n1 : test_ind) {
	    	 	   std::cout << n1 << ' ';
	    		}
				std::cout << "\n" << std::flush;
				/*
				file_performance << "\ntrain: \t";
				for(int n1 : train_ind) {
	    	 	   file_performance << n1 << ' ';
	    		}
				file_performance << "\ntest: \t";
				for(int n1 : test_ind) {
	    	 	   file_performance << n1 << ' ';
	    		}
				file_performance << "\n" << std::flush;
				*/
				
				for(int n1 : train_ind) {
	    	 	   file_performance << n1 << ", ";
	    		}
				
				
				
				
				
				
				
				// Least Sq:
				// Change 
				int do_least_sq = 1;	// 0 Subrata -- least sq have better initial likelihood-but stucks and gives nan in 	some value
				Matrix_eig_row W_init = Init_val(train, TE_train, TR_train, our_dim_train, 
				                             r_scale, TE_scale, TR_scale, W1_init, W2_init, do_least_sq, will_write);
				Debug1("W_init after LS: ");
				show_head(W_init);
				std::cout << std::flush;
				
				
				perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
				perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
				perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
				perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
				Debug0("Avg perfs LS: " << perf_1.mean() << ", " << perf_2.mean() << ", "
									 << perf_3.mean() << ", " << perf_4.mean());
				std::cout << "Performances over images LS: " << perf_1.transpose() << "\n";
				std::cout << "Performances over images LS: " << perf_2.transpose() << "\n";
				std::cout << "Performances over images LS: " << perf_3.transpose() << "\n";
				std::cout << "Performances over images LS: " << perf_4.transpose() << "\n" << std::endl;

				file_performance << perf_1.mean() << ", " << perf_2.mean() << ", "
									 << perf_3.mean() << ", " << perf_4.mean() << ", ";
//				file_performance << perf_1.mean() << ", " << perf_2.mean() << ", "
//									 << perf_3.mean() << ", " << perf_4.mean() << "\n";
//				file_performance << "Avg perfs LS: " << perf_1.mean() << ", " << perf_2.mean() << ", "
//									 << perf_3.mean() << ", " << perf_4.mean() << "\n"; 
//				file_performance << "Performances over images LS: \t" << perf_1.transpose() << "\n";
//				file_performance << "Performances over images LS: \t" << perf_2.transpose() << "\n";
//				file_performance << "Performances over images LS: \t" << perf_3.transpose() << "\n";
//				file_performance << "Performances over images LS: \t" << perf_4.transpose() << "\n\n\n";
				
				
				
			
				
				MRF_param MRF_obj_1(our_dim_train[1], our_dim_train[2], our_dim_train[3]);
				
				
				
				Matrix_eig_row W_LS = W_init;
				
				
				
				// Likelihood Based optimization:
				
				Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
				Vector_eig beta_init = 1.0*Vector_eig::Ones(3);
			
					
				
				// Non -penalized:
				
				OSL_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
				          our_dim_train[1], our_dim_train[2], our_dim_train[3], r_scale, TE_scale, TR_scale, MRF_obj_1, 
				          100, 0, 0.1, 1e-5, 1);
				//change
				show_head(W_init);
				
				
				
				
				
				perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
				perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
				perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
				perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
				Debug0("Avg perfs MLE: " << perf_1.mean() << ", " << perf_2.mean() << ", "
									 << perf_3.mean() << ", " << perf_4.mean());
				std::cout << "Performances over images Likelihood: " << perf_1.transpose() << "\n";
				std::cout << "Performances over images Likelihood: " << perf_2.transpose() << "\n";
				std::cout << "Performances over images Likelihood: " << perf_3.transpose() << "\n";
				std::cout << "Performances over images Likelihood: " << perf_4.transpose() << "\n" << std::endl;
				
				file_performance << perf_1.mean() << ", " << perf_2.mean() << ", "
									 << perf_3.mean() << ", " << perf_4.mean() << ", ";
//				file_performance << perf_1.mean() << ", " << perf_2.mean() << ", "
//									 << perf_3.mean() << ", " << perf_4.mean() << "\n";
//				file_performance << "Avg perfs MLE: " << perf_1.mean() << ", " << perf_2.mean() << ", "
//									 << perf_3.mean() << ", " << perf_4.mean() << "\n";
//				file_performance << "Performances over images Likelihood: \t" << perf_1.transpose() << "\n";
//				file_performance << "Performances over images Likelihood: \t" << perf_2.transpose() << "\n";
//				file_performance << "Performances over images Likelihood: \t" << perf_3.transpose() << "\n";
//				file_performance << "Performances over images Likelihood: \t" << perf_4.transpose() << "\n\n\n";
				
				
				
				
				W_init.noalias() = W_LS;
				
				
				
				
				// Penalised:
				
				OSL_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
				          our_dim_train[1], our_dim_train[2], our_dim_train[3], r_scale, TE_scale, TR_scale, MRF_obj_1, 
				          100, 1, 0.1, 1e-5, 1);
				//change
				Debug1("W - Penalized Likelihood");
				show_head(W_init);
				
				
				
				
				
				perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
				perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
				perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
				perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
				Debug0("Avg perfs MPLE: " << perf_1.mean() << ", " << perf_2.mean() << ", "
									 << perf_3.mean() << ", " << perf_4.mean());
				std::cout << "Performances over images Penalized: " << perf_1.transpose() << "\n";
				std::cout << "Performances over images Penalized: " << perf_2.transpose() << "\n";
				std::cout << "Performances over images Penalized: " << perf_3.transpose() << "\n";
				std::cout << "Performances over images Penalized: " << perf_4.transpose() << "\n" << std::endl;
				
				file_performance << perf_1.mean() << ", " << perf_2.mean() << ", "
									 << perf_3.mean() << ", " << perf_4.mean() << "\n";
				std::cout << std::flush;
//				file_performance << perf_1.mean() << ", " << perf_2.mean() << ", "
//									 << perf_3.mean() << ", " << perf_4.mean() << "\n";
//				file_performance << "Avg perfs MPLE: " << perf_1.mean() << ", " << perf_2.mean() << ", "
//									 << perf_3.mean() << ", " << perf_4.mean() << "\n";
//				file_performance << "Performances over images Penalized: \t" << perf_1.transpose() << "\n";
//				file_performance << "Performances over images Penalized: \t" << perf_2.transpose() << "\n";
//				file_performance << "Performances over images Penalized: \t" << perf_3.transpose() << "\n";
//				file_performance << "Performances over images Penalized: \t" << perf_4.transpose() << "\n\n\n";
			}
			
		}
//	}
	
	
	
	
	
	
	
	
	
	
	
	
	file_performance.close();
	
	
	
	
	
	
	std::time_t t2 = std::time(nullptr);
	std::tm tm2 = *std::localtime(&t2);
	std::cout << "Current time: " << std::put_time(&tm2, "%c %Z") << '\n';

	return 0;
}

