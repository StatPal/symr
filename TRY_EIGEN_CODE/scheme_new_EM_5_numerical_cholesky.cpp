/**
*
* See this:
g++ scheme_new_EM_5_numerical_cholesky.cpp -o test -I /usr/include/eigen3 -O3
g++ ~/MRI/Headers/TRY_EIGEN_2/scheme_new_EM_5_numerical.cpp -o scheme_new_EM_5_numerical -I ~/MRI/Headers -O3 -std=c++11
*
./test ../Read_Data/ZHRTS1.nii Dummy_sd.txt 0
./scheme_new_EM_5_numerical ../data/ZHRTS1.nii ../data/Dummy_sd.txt 0 > scheme_new_EM_5_numerical.txt


./test ../Read_Data/new_phantom.nii Dummy_sd.txt 0






* Notes:
Subrata -- do Cholesky -- done
In gradient also!!  -- done


There is negative in l_star() - But not in grad vector
Subrata - correct it. -- done

Check sign of grad - is there opposite sign in MRF part ?
Because, Psi_inv(0,0) increases with iterations - and objective function worsens?
Checked


*/



#include "scheme_new_numerical.hpp"
#include "Read_files_2.hpp"
#include "Init_value_6_numerical.hpp"

#include "../CppNumericalSolvers/include/cppoptlib/meta.h"
#include "../CppNumericalSolvers/include/cppoptlib/boundedproblem.h"
#include "../CppNumericalSolvers/include/cppoptlib/solver/lbfgsbsolver.h"





const Matrix_eig G((Matrix_eig(6,9) << 
  1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1).finished());


/***************************************************
************** Function Writing part ***************
****************************************************/


/*
Penalised Negative log likelihood -- to be minimised
Matrix sizes: nx3, 3x3, 3(2)x1, mx1, mx1, mx1, nxm, ...
No change for cholesky inside the function
*/
double l_star(const Matrix_eig &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta,
              const Vector_eig &TE, const Vector_eig &TR, 
              const Vector_eig &sigma, const Matrix_eig &r, 
              int n_x, int n_y, int n_z){

	Matrix_eig v = v_mat(W, TE, TR);
	int m = v.cols();
	int n = v.rows();
	double tmp2 = 0.0, tmp3 = 0.0, tmp1 = 0.0;

	//Rice part://
	long double likeli_sum = 0.0, tmp4 = 0.0;
	auto start = std::chrono::high_resolution_clock::now();
	int i = 0, j = 0;
	for(i = 0; i < n; ++i) {
		for(j = 0; j < m; ++j) {
			tmp2 = r(i,j)/SQ(sigma(j));
			tmp3 = (SQ(r(i,j))+SQ(v(i,j)))/SQ(sigma(j));
			tmp1 = logBesselI0(tmp2*v(i,j));
			likeli_sum += (log(tmp2) + tmp1 - 0.5*tmp3) ;
		}
	}
	
	//MRF part://
	double tmp = (Psi_inv*W.transpose()*Lambda(beta, n_x, n_y, n_z)*W).trace();
	likeli_sum += ( -tmp + 3*sp_log_det_specific(beta, n_x, n_y, n_z) + n*log_det_3(Psi_inv) - 3*n*log(2*M_PI) )/2;
	
	return -likeli_sum;
}


/*
* Negative Penalised Q function - to be minimised
* Matrix sizes: nx3, 3x3, 3(2)x1, mx1, mx1, mx1, nxm, ...
* Checks are removed as much as possible to get a cleaner code
* Add them whenever necessary (see old git commits)
*/
double Q_star(const Matrix_eig &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta,
              const Vector_eig &TE, const Vector_eig &TR, 
              const Vector_eig &sigma, const Matrix_eig &r, const Matrix_eig &W_old,
              int n_x, int n_y, int n_z){

	
	Matrix_eig v = v_mat(W, TE, TR);
	Matrix_eig v_old = v_mat(W_old, TE, TR); 
	int m = v.cols();
	int n = v.rows();
	
	double tmp2 = 0.0, tmp3 = 0.0;
	
	//Rice part://
	double likeli_sum = 0.0;
	int i = 0, j = 0;
	for(i = 0; i < n; ++i) {
		for(j = 0; j < m; ++j) {
			tmp2 = r(i,j)*v_old(i,j)/SQ(sigma(j));
			tmp3 = besselI1_I0(tmp2);
			likeli_sum += v(i,j)*(- 0.5*v(i,j) + r(i,j)*tmp3)/SQ(sigma(j));
		}
	}
	
	//if(std::isnan(likeli_sum)){
	//	Debug0("nan after Rice part");
	//}
	
	// MRF part://
	double tmp = (Psi_inv*W.transpose()*Lambda(beta, n_x, n_y, n_z)*W).trace();
	likeli_sum += ( -tmp + 3*sp_log_det_specific(beta, n_x, n_y, n_z) + n*log_det_3(Psi_inv) - 3*n*log(2*M_PI) )/2;
	//if(std::isnan(likeli_sum))
	//	Debug0("nan after MRF part");
	
	Debug2(" - Q function calculated: " << -likeli_sum);
	
	return (-likeli_sum);
}






/*
* Negative Gradient of Penalised Q function
*/
Vector_eig Q_grad_vec(const Matrix_eig &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                      const Vector_eig &TE, const Vector_eig &TR, 
                      const Vector_eig &sigma, const Matrix_eig &r, const Matrix_eig &W_old,
                      int n_x, int n_y, int n_z){

	
	Matrix_eig v = v_mat(W, TE, TR);
	Matrix_eig v_old = v_mat(W_old, TE, TR); 
	int n = n_x * n_y * n_z;
	int m = v.cols();
	double temp = 0.0, tmp2 = 0.0, tmp3 = 0.0;
	SpMat Gamma_inv = Lambda(beta, n_x, n_y, n_z);
	Matrix_eig W_grad(n, 3);
	
	// MRF contribution part//
	Matrix_eig MRF_grad = Gamma_inv * W * Psi_inv;
	
	// W - grad part//
	int i = 0, j = 0, k = 0;
	for(i = 0; i < n; ++i) {
		for(k = 0; k < 3; ++k) {
			temp = 0.;
			for(j = 0; j < m ; ++j){
				tmp2 = r(i,j)/SQ(sigma(j));
				tmp3 = -v(i,j)/SQ(sigma(j)) + tmp2*besselI1_I0(tmp2*v_old(i,j));
				temp += tmp3* simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k);
			}
			W_grad(i, k) = temp - MRF_grad(i,k);
		}
	}
	
	
	
	// Other - grad part//
	Matrix_eig Psi_grad = 0.5 * G * to_vector(n * Psi_inv.llt().solve(Matrix3d_eig::Identity(3, 3)) - W.transpose()*Gamma_inv*W);
	
	Matrix_eig temp_mat = W * Psi_inv;
	double beta_x_grad = 1.5*sp_log_inv_specific(beta, n_x, n_y, n_z, 0) - 
	                     0.5*(W.transpose()*Kron_Sparse_eig(J_n(n_x), I_n(n_y*n_z))*temp_mat).trace();
	double beta_y_grad = 1.5*sp_log_inv_specific(beta, n_x, n_y, n_z, 1) - 
	                     0.5*(W.transpose()*Kron_Sparse_eig(Kron_Sparse_eig(I_n(n_x), J_n(n_y)), I_n(n_z))*temp_mat).trace();
	
	
	Debug2("Grad calculated.\n");
	return(-to_param_vec_grad(W_grad, Psi_grad, beta_x_grad, beta_y_grad));	// Check this fn and EM_14
}
// -ve sign? --Subrata
// No change for cholesky inside the function - change in optimizer















/***************************************************
**************** Optimization part *****************
****************************************************/

/*
This is a minimiser!
*/
template<typename T>
class EM_opt : public cppoptlib::BoundedProblem<T> {
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	TMatrix r;
	TVector r2;


  public:
	EM_opt(const TVector y_) :											// Is it constructor?
		cppoptlib::BoundedProblem<T>(y_.size()), r2(y_){}




	int n_x, n_y, n_z;
	double beta_z = 0.1;												//Subrata - or get the value. 
	TVector TE, TR, sigma;
	TMatrix W_old;
	
	
	// Track the best:
	TVector current_best_param;
	T current_best_val = 1.0e+15;
	
	
	// Extra check:
	TVector lb, ub;
	double prev_val = 0.0, prev_lik = 0.0;								// necessary??




	
	// all_param(8) has the following format:
	// First 3*n elements for W! (0 to n-1 is first column of W etc.)
	// Then 6 elements for Cholesky decomposition of Psi_inv!
	// Then 2 element for beta_x and beta_y!
	
	T value(const TVector &all_param) {
	
		
		int n = n_x*n_y*n_z;
		
		//W
		Matrix_eig W = Matrix_eig::Zero(n, 3);
		W.col(0) = all_param.segment(  0, n);
		W.col(1) = all_param.segment(  n, n);
		W.col(2) = all_param.segment(2*n,n);
		check_bounds(W, lb, ub);			// I guess not working - take a look later
		
		//show_head(W, 2);					// What is happening!
		//show_head_vec(all_param, 10, 1);
		
		
		//Psi
		Vector_eig temp_L = all_param.segment(3*n, 6);
		Matrix3d_eig L_mat = to_L_mat(temp_L); 				// 6 element vec to Lower traingular mat: [a0, a1, a2; 0, a3, a4; 0, 0, a5]
		Matrix3d_eig Psi_inv = from_Cholesky(L_mat);		// Psi from Cholesky decomposition!
		
		//beta
		Vector_eig beta = Vector_eig::Zero(3);
		beta(0) = all_param(3*n+6); beta(1) = all_param(3*n+7); beta(2) = beta_z;
		
		
		
		double temp_val = Q_star(W, Psi_inv, beta, TE, TR, sigma, r, W_old, n_x, n_y, n_z);				// Bug - is it still there? Hey bug!!
		double temp_lik = l_star(W, Psi_inv, beta, TE, TR, sigma, r, n_x, n_y, n_z);
		Debug1("- Q fn:" << temp_val);
		Debug1("- log Likelihood: " << temp_lik);
		
		
		// Track immediate past:
		double diff_val = prev_val - temp_val;			// smaller temp_val is better. 
		double diff_lik = prev_lik - temp_lik;
		
		if(diff_val > 0.0){
			Debug0("Better than previous value. Previous val:"<< prev_val << 
					"; current_val:" << temp_val << "; diff:" << diff_val);
		} else if(diff_val == 0){
			Debug1("Same Q");
		} else if( diff_val < 0.0){
			Debug1("Increase, bad, difference: " << diff_val);
		}
		Debug1("Diff in likelihood: "<< diff_lik << "\n");
		Debug2("Previous val:"<< prev_val << "; current_val:" << temp_val);
		
		prev_val = temp_val;
		prev_lik = temp_lik;
		
		
		
		// Track the best (among all the loops, not just the past):
		if(temp_val < current_best_val){
			current_best_param = all_param;
			current_best_val = temp_val;
		}
		
		return (temp_val);
	}




	
	// all_param(8) has the following format:
	// First 3*n elements for W! (0 to n-1 is first column of W etc.)
	// Then 6 elements for Cholesky decomposition of Psi_inv!
	// Then 2 element for beta_x and beta_y!
	
	void gradient(const TVector &all_param, TVector &Y) {
	
		
		//check_nan_vec(all_param);
		int n = n_x*n_y*n_z;
		
		//W
		Matrix_eig W = Matrix_eig::Zero(n, 3);
		W.col(0) = all_param.segment(0,  n);
		W.col(1) = all_param.segment(n,  n);
		W.col(2) = all_param.segment(2*n,n);
		//check_bounds(W, lb, ub);
		
		//Psi
		Vector_eig temp_L = all_param.segment(3*n, 6);
		Matrix3d_eig L_mat = to_L_mat(temp_L);
		Matrix3d_eig Psi_inv = from_Cholesky(L_mat);
		
		//beta
		Vector_eig beta = Vector_eig::Zero(3);
		beta(0) = all_param(3*n+6); beta(1) = all_param(3*n+7); beta(2) = beta_z;
		
		
		Y = Q_grad_vec(W, Psi_inv, beta, TE, TR, sigma, r, W_old, n_x, n_y, n_z);
		
		
		
		// Chain rule adjustment for Cholesky:
		Vector_eig chain = Y.segment(3*n, 6);
		Y.segment(3*n, 6) = to_grad_Cholesky(temp_L)*chain;		//check transpose
		
		//Debug3("grad: " << Y.transpose() << "\n" );
		
	}

};







/*
* Main optimizaing fn for direct EM.
*/
void EM_solve(Matrix_eig &W_init, Matrix3d_eig &Psi_inv, Vector_eig &beta, 
              const Vector_eig &TE_example, const Vector_eig &TR_example, 
              const Vector_eig &sigma, const Matrix_eig &r,
              int n_x, int n_y, int n_z, double TE_scale, double TR_scale, 
              int maxiter = 20, double abs_diff = 1e-6, int verbose = 0){




	
	Debug0("Doing EM Estimate!");
	auto time_1_likeli = std::chrono::high_resolution_clock::now();
	
	int iter = 0;
	int n = n_x * n_y * n_z;
	EM_opt<double> f(Vector_eig::Ones(3));
	





	
	// Initial other parameters (e.g., dimension): //
	
	
	f.n_x = n_x; f.n_y = n_y; f.n_z = n_z;
	f.sigma = sigma;	f.r = r;	f.TE = TE_example;	f.TR = TR_example;
	Matrix_eig W_old = W_init;
	f.W_old = W_old;
	
	
	
	
	// Bounds: //
	
	Vector_eig lb = Vector_eig::Ones(3*n+6+1+1, 1);
	Vector_eig ub = Vector_eig::Ones(3*n+6+1+1, 1);
	
	lb.segment(  0, n) *= 0.0;
	lb.segment(  n, n) *= exp(-1/(0.01*TR_scale));
	lb.segment(2*n, n) *= exp(-1/(0.001*TE_scale));
	ub.segment(  0, n) *= 450;
	ub.segment(  n, n) *= exp(-1/(4.0*TR_scale));
	ub.segment(2*n, n) *= exp(-1/(0.2*TE_scale));
	
	
	// Bound on Cholesky etc: (Why 255 btw? - Forgot)
	ub.segment(3*n, 6+1+1) *= 255;
	lb.segment(3*n, 6+1+1) = -255*Vector_eig::Ones(6+1+1);
	lb(3*n) = 1e-5; lb(3*n+3) = 1e-5;lb(3*n+5) = 1e-5;	lb(3*n+6) = 1e-5;lb(3*n+7) = 1e-5;
	//ub.segment(3*n, 6+1+1) *= 15;		// Get's sqaured
	//lb.segment(3*n, 6+1+1) = 0.00001*Vector_eig::Ones(6);
	
	f.setLowerBound(lb);	f.setUpperBound(ub);
	f.ub = ub; 		f.lb = lb;
	
	







	
	// Optimization parameter:
	
	Vector_eig param_new = Vector_eig::Ones(3*n+6+1+1, 1);
	
	//Psi_inv = Eigen::Matrix3d::Identity()*1.1;	// 1.1 is just for test of ub
	Matrix3d_eig L( Psi_inv.llt().matrixL() );
	param_new.segment(3*n, 6) = from_L_mat(L);
	param_new(3*n + 6) = beta(0); param_new(3*n + 7) = beta(1);
	Vector_eig param_old = param_new;
	
	// Vector_eig param_new = to_param_vec(W_init, Psi_inv, beta_x, beta_y);	//W_init.row(i); // BUG in cholesky part  - confused
	Debug0 ("f(x) at first: \n");
	Debug0(f.value(param_new));					// Bug? -- are you still there?
	







	
	double old_val = 1.0e+15, old_likeli = 1.0e+15, current_best_likeli = 1.0e+15;
	int bad_count_o = 0, bad_count_o_2 = 0, bad_bound_1 = 0, bad_bound_2 = 0;
	
	
	// Declaring the solver:
	cppoptlib::LbfgsbSolver<EM_opt<double>> solver;
	
	// Subrata - Setting the parameters: new  -- (see simple_withoptions.cpp)
	cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults(); 	// Create Criteria class to set the stopping rule
	//crit.iterations = 10000;														// Change the number of allowed iterations
	//solver.setStopCriteria(crit);
	
	
	
	//* EM loop *//
	while(iter < maxiter){
	
	
		Debug1("\n--------------------------------------------------------------------------\nIteration: " << iter++ << "\n");
		
		W_old.col(0) = param_old.segment(0,  n);
		W_old.col(1) = param_old.segment(n,  n);
		W_old.col(2) = param_old.segment(2*n,n);
		f.W_old = W_old;
		old_val = f.value(param_old);
		
		
		// Track the best:
		f.current_best_val = 1.0e+15;
		
		//Solve:
		solver.minimize(f, param_new);
		Debug1("Solved for this iteration!\n");
		
		// Track the best:
		param_new = f.current_best_param;
		double best_val = f.current_best_val;
		
		Debug0("f(param_new) in argmin: " << best_val << "\t old val:" << old_val);
		std::cout << "Solver status: " << solver.status() << "\n";	//Guess: bad reports: under constraints => grad is not ~0 
		std::cout << "Final criteria values: " << "\n" << solver.criteria() << "\n";
		
		
		//moving on part:
		if(abs_sum(param_old - param_new) <= abs_diff){
			std::cout << "Stopped after " << iter << " iterations" << "\n";
			break;
		}
		
		if(best_val >= old_val){ 
			Debug1("Value not decreased!! old val: " << old_val << ";\t new val: " << best_val);
			bad_count_o++;
		}
		if(best_val > old_val){
			bad_count_o_2++;
		}
		//if(best_val<old_val){		// Good case? 
		// -- without this, diff between abs sum between old and new parameters would be 0!
			param_old = param_new;
		//}
		
		Debug1("Another iteration done\n\n\n");
	}
	
	
	
	//W
	W_init.col(0) = param_new.segment(0, n);
	W_init.col(1) = param_new.segment(n, n);
	W_init.col(2) = param_new.segment(2*n,n);
	check_bounds(W_init, lb, ub);
	
	//Psi
	Vector_eig temp_L = param_new.segment(3*n, 6);
	Matrix3d_eig L_mat = to_L_mat(temp_L);
	Psi_inv = from_Cholesky(L_mat);
	
	//beta
	beta(0) = param_new(3*n+6); beta(1) = param_new(3*n+7); beta(2) = f.beta_z;
	
	
	
	auto time_end_likeli = std::chrono::high_resolution_clock::now();
	auto duration_likeli = std::chrono::duration_cast<std::chrono::seconds>(time_end_likeli - time_1_likeli);
	Debug1("Time taken for whole EM: " << duration_likeli.count() << " seconds\n");
	
}











int main(int argc, char * argv[]) {
	
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
	double TE_scale = 1.01/TE_example.minCoeff();
	double TR_scale = 1.01/TR_example.minCoeff();
	Debug0("TE scale: " << TE_scale);
	Debug0("TR scale: " << TR_scale);
	TE_example *= TE_scale;
	TR_example *= TR_scale;
	
	/*
	Vector_eig lb(3), ub(3);
	lb << 0.0, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	*/
	
	double W1_init = exp(-1/(2.0*TR_scale));		// exp(-1/(2.0*1.01))
	double W2_init = exp(-1/(0.1*TE_scale));		// exp(-1/(0.1*1.01/0.03))
	





	
	Matrix_eig r = Preprocess_data(data_file, our_dim, will_write);
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
	for(int i = 0; i < train_ind.size(); ++i){
		train.col(i) = r.col(train_ind[i]);
		TE_train[i] = TE_example(train_ind[i]);
		TR_train[i] = TR_example(train_ind[i]);
		sigma_train[i] = sigma(train_ind[i]);
	}
	for(int i = 0; i < 8; ++i){
		our_dim_train[i] = our_dim[i];
	}
	our_dim_train[4] = (short)train_ind.size();		//our_dim[0] = 3 or 4
	
	Matrix_eig test(r.rows(), test_ind.size());
	Vector_eig TE_test(test_ind.size()), TR_test(test_ind.size()), sigma_test(test_ind.size());
	for(int i = 0; i < test_ind.size(); ++i){
		test.col(i) = r.col(test_ind[i]);
		TE_test[i] = TE_example(test_ind[i]);
		TR_test[i] = TR_example(test_ind[i]);
		sigma_test[i] = sigma(test_ind[i]);
	}
	
	
	
	// Temp result: Performance on the Init W: 
	Matrix_eig W_1st = Matrix_eig::Ones(our_dim_train[1]*our_dim_train[2]*our_dim_train[3], 3);
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
		if(W_1st(i, 0)>450){
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
	// change
	int do_least_sq = 1;	// 0 Subrata -- least sq have better initial likelihood-but stucks and gives nan in some value
	Matrix_eig W_init = Init_val(train, TE_train, TR_train, our_dim_train, 
	                             TE_scale, TR_scale, W1_init, W2_init, do_least_sq, will_write);
	Debug2("W initial done");
	check_nan(W_init);
	show_head(W_init);
	perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
	std::cout << "Performances over images: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images: " << perf_4.transpose() << "\n";
	





	
	// Likelihood Based optimization:
	int n = our_dim[1]*our_dim[2]*our_dim[3];
	Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
	Vector_eig beta_init = 0.1*Vector_eig::Ones(3);
	EM_solve(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma, train, 
	         our_dim_train[1], our_dim_train[2], our_dim_train[3], TE_scale, TR_scale);
	show_head(W_init);
	
	
	perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
	std::cout << "Performances over images: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images: " << perf_4.transpose() << "\n";
	

	//Hessian_mat_without_MRF(W_init, Psi_inv_init, beta_init, TE_example,TR_example, sigma, r, our_dim[1], our_dim[2], our_dim[3]);

	return 0;
}


