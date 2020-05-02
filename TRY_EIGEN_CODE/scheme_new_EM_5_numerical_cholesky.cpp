/*
*
* See this:
g++ scheme_new_EM_5_numerical_cholesky.cpp -o test -I /usr/include/eigen3 -O3 --std=c++17
g++ ~/MRI/Headers/TRY_EIGEN_2/scheme_new_EM_5_numerical.cpp -o scheme_new_EM_5_numerical -I ~/MRI/Headers -O3 -std=c++11
* 
*
./test ../Read_Data/ZHRTS1.nii Dummy_sd.txt 0
./test ../data/ZHRTS1.nii ../data/Dummy_sd.txt 0
./scheme_new_EM_5_numerical ../data/ZHRTS1.nii ../data/Dummy_sd.txt 0 > scheme_new_EM_5_numerical.txt
*/



#include "scheme_new_numerical.hpp"
#include "Read_files_2.hpp"
#include "Init_value_6_numerical.hpp"

#include "../optim_cpp_solver/include/cppoptlib/meta.h"
#include "../optim_cpp_solver/include/cppoptlib/boundedproblem.h"
#include "../optim_cpp_solver/include/cppoptlib/solver/lbfgsbsolver.h"


//Subrata 
// Make it --std=c++11 -- remove bessel std::cyl_bessel_i - do later


//Subrata -- do Cholesky 
// -- in all proper places -- 
// In gradient also!!


// There is negative in l_star() - But not in grad vector
// Subrata - correct it. 

// Check sign of grad - is there opposite sign in MRF part ?
// Because, Psi_inv(0,0) increases with iterations - and objective function worsens


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;

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
*/
double l_star(Matrix_eig W, Matrix3f_eig Psi_inv, Vector_eig beta,
              Vector_eig TE, Vector_eig TR, Vector_eig sigma, Matrix_eig r, 
              int n_x, int n_y, int n_z){	// nx3, 3 or 2, mx1...

	Matrix_eig v = v_mat(W, TE, TR);
	int m = v.cols();
	int n = v.rows();
	double tmp2 = 0.0, tmp3 = 0.0, tmp1 = 0.0;

	//Rice part://
	long double likeli_sum = 0.0, tmp4 = 0.0;
	auto start = std::chrono::high_resolution_clock::now();
	for(long int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
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
// No change for cholesky inside the function

/*
* Negative Penalised Q function - to be minimised
*/
double Q_star(Matrix_eig W, Matrix3f_eig Psi_inv, Vector_eig beta,
              Vector_eig TE, Vector_eig TR, Vector_eig sigma, Matrix_eig r, Matrix_eig W_old,
              int n_x, int n_y, int n_z){	// nx3, 3 or 2, mx1...
	
	Debug2("Calculation of Q starts:");
	check_nan(W);
	Matrix_eig v = v_mat(W, TE, TR);
	Matrix_eig v_old = v_mat(W_old, TE, TR); 
	int m = v.cols();
	int n = v.rows();
	Debug2("Values of m, n: " << m << ", " << n);
	check_nan(v);
	check_nan(v_old);
	check_nan(r);
	
	double tmp2 = 0.0, tmp3 = 0.0;
	
	/*Rice part:*/
	double likeli_sum = 0.0;
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
			tmp2 = r(i,j)*v_old(i,j)/SQ(sigma(j));
			//if(std::isnan(tmp2)){Debug1("r(i,j): " << r(i,j) << " v_old:" << v_old(i,j) << " sigma^2: " << SQ(sigma(j)));}
			tmp3 = besselI1_I0(tmp2);				// Mistake --Subrata  --corrected
			//if(std::isnan(tmp3)){Debug1("tmp3:" << tmp3);}
			likeli_sum += v(i,j)*(- 0.5*v(i,j) + r(i,j)*tmp3)/SQ(sigma(j));
		}
	}
	
	if(std::isnan(likeli_sum)){
		Debug0("nan after Rice part");
	}
	
	/*MRF part:*/
	double tmp = (Psi_inv*W.transpose()*Lambda(beta, n_x, n_y, n_z)*W).trace();
	if(std::isnan(tmp)){
		Debug0("nan at tmp part");
	}
	likeli_sum += ( -tmp + 3*sp_log_det_specific(beta, n_x, n_y, n_z) + n*log_det_3(Psi_inv) - 3*n*log(2*M_PI) )/2;
	
	if(std::isnan(log_det_3(Psi_inv))){
		Debug0("nan at log_det_3: Psi_inv:\n"<< Psi_inv << "\nDet:" << log_det_3(Psi_inv));
	}
	Debug1("Psi_inv: \n" << Psi_inv);
	Debug1("beta: " << beta.transpose());
	
	if(std::isnan(likeli_sum)){
		Debug0("nan after MRF part");
	}
	
	check_nan(W);
	Debug2(" - Q function calculated: " << -likeli_sum);
	
	
	assert( ! std::isnan(-likeli_sum) );		//Better way?
	
	return (-likeli_sum);
}
// -ve sign? -- Subrata
// No change for cholesky inside the function






/*
* Negative Gradient of Penalised Q function
*/
Vector_eig Q_grad_vec(Matrix_eig W, Matrix3f_eig Psi_inv, Vector_eig beta, 
                   Vector_eig TE, Vector_eig TR, Vector_eig sigma, Matrix_eig r, Matrix_eig W_old,
                   int n_x, int n_y, int n_z){
	
	Debug2("Grad calculation started!");
	check_nan(W);
	Matrix_eig v = v_mat(W, TE, TR);
	Matrix_eig v_old = v_mat(W_old, TE, TR); 
	int n = n_x * n_y * n_z;
	int m = v.cols();
	double temp = 0.0, tmp2 = 0.0, tmp3 = 0.0;
	SpMat Gamma_inv = Lambda(beta, n_x, n_y, n_z);
	Matrix_eig W_grad = W;
	
	/* MRF contribution part*/
	Matrix_eig MRF_grad = Gamma_inv * W * Psi_inv;
	
	/* W - grad part*/
	for(int i = 0; i < n; ++i){
		for(int k = 0; k < 3; ++k){
			temp = 0.;									// Missed this -- correct this.
			for(int j = 0; j < m ; ++j){
				tmp2 = r(i,j)/SQ(sigma(j));
				tmp3 = -v(i,j)/SQ(sigma(j)) + tmp2*besselI1_I0(tmp2*v_old(i,j));
				temp += tmp3* simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k);
			}
			W_grad(i, k) = temp - MRF_grad(i,k);
		}
	}
	Debug2("Grad rice part calculated.\n");
	
	/* Other - grad part*/
	Matrix_eig Psi_grad = 0.5 * G * to_vector(n * Psi_inv.llt().solve(Matrix3f_eig::Identity(3, 3)) - W.transpose()*Gamma_inv*W);
	
	Matrix_eig temp_mat = W * Psi_inv;
	double beta_x_grad = 1.5*sp_log_inv_specific(beta, n_x, n_y, n_z, 0) - 0.5*(W.transpose()*Kron_Sparse_eig(J_n(n_x), I_n(n_y*n_z))*temp_mat).trace();
	double beta_y_grad = 1.5*sp_log_inv_specific(beta, n_x, n_y, n_z, 1) - 0.5*(W.transpose()*Kron_Sparse_eig(Kron_Sparse_eig(I_n(n_x), J_n(n_y)), I_n(n_z))*temp_mat).trace();
	
	check_nan(W);
	Debug2("Grad calculated.\n");
	return(-to_param_vec_grad(W_grad, Psi_grad, beta_x_grad, beta_y_grad));
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

	TVector TE, TR, sigma;
	int n_x, n_y, n_z;
	double beta_z = 0.1;												//Subrata - or get the value. 
	TMatrix W_old;
	
	// Track the best:
	TVector current_best_param;
	T current_best_val = 1.0e+15;
	
	// Extra check:
	TVector lb, ub;
	double prev_val = 0.0, prev_lik = 0.0;


	T value(const TVector &all_param) {
	
		Debug2("Inside computation of value: ");
		check_nan_vec(all_param);
		int n = n_x*n_y*n_z;
		
		//W
		Matrix_eig W = Matrix_eig::Zero(n, 3);
		W.col(0) = all_param.segment(0, n);
		W.col(1) = all_param.segment(n, n);
		W.col(2) = all_param.segment(2*n,n);
		check_bounds(W, lb, ub);
		
		//Psi
		Vector_eig temp_L = all_param.segment(3*n, 6);
		Matrix3f_eig L_mat = to_L_mat(temp_L);
		Matrix3f_eig Psi_inv = from_Cholesky(L_mat);
		
		//beta
		Vector_eig beta = Vector_eig::Zero(3);
		beta(0) = all_param(3*n+6); beta(1) = all_param(3*n+7); beta(2) = beta_z;
		
		double temp_val = Q_star(W, Psi_inv, beta, TE, TR, sigma, r, W_old, n_x, n_y, n_z);
		double temp_lik = l_star(W, Psi_inv, beta, TE, TR, sigma, r, n_x, n_y, n_z);
		Debug1("- Q fn:" << temp_val);
		Debug1("- log Likelihood: " << temp_lik);
		
		
		// Track immediate past:
		double diff_val = prev_val - temp_val;			// smaller temp_val is better. 
		double diff_lik = prev_lik - temp_lik;
		
		if(diff_val > 0.0){
			Debug0("Better than previous value. Previous val:"<< prev_val << "; current_val:" << temp_val << "; diff:" << diff_val);
		} else if(diff_val == 0){
			Debug1("Same Q");
		} else if( diff_val < 0.0){
			Debug1("Increase, bad, difference: " << diff_val);
		}
		Debug1("Diff in likelihood: "<< diff_lik);
		Debug2("Previous val:"<< prev_val << "; current_val:" << temp_val);
		
		prev_val = temp_val;
		prev_lik = temp_lik;
		
		
		
		// Track the best:
		if(temp_val < current_best_val){
			current_best_param = all_param;
			current_best_val = temp_val;
		}
		
		Debug1("Computation of values done!");
		return (temp_val);
	}

	void gradient(const TVector &all_param, TVector &Y) {
	
		Debug2("Inside gradient: ");
		check_nan_vec(all_param);
		int n = n_x*n_y*n_z;
		
		//W
		Matrix_eig W = Matrix_eig::Zero(n, 3);
		W.col(0) = all_param.segment(0,  n);
		W.col(1) = all_param.segment(n,  n);
		W.col(2) = all_param.segment(2*n,n);
		check_bounds(W, lb, ub);
		
		//Psi
		Vector_eig temp_L = all_param.segment(3*n, 6);
		Matrix3f_eig L_mat = to_L_mat(temp_L);
		Matrix3f_eig Psi_inv = from_Cholesky(L_mat);
		
		//beta
		Vector_eig beta = Vector_eig::Zero(3);
		beta(0) = all_param(3*n+6); beta(1) = all_param(3*n+7); beta(2) = beta_z;
		
		Y = Q_grad_vec(W, Psi_inv, beta, TE, TR, sigma, r, W_old, n_x, n_y, n_z);
		
		// Chain rule for Cholesky:
		Vector_eig chain = Y.segment(3*n, 6);
		Y.segment(3*n, 6) = to_grad_Cholesky(temp_L)*chain;		//check transpose
		
		check_nan_vec(Y);
		Debug1("Computation of grad done!\n");
		Debug2("grad: " << Y.transpose() << "\n" );
	}
};


void EM_solve(Matrix_eig &W_init, Matrix3f_eig &Psi_inv, Vector_eig &beta, 
                Vector_eig TE_example, Vector_eig TR_example, Matrix_eig r, 
                int n_x, int n_y, int n_z, Vector_eig sigma, int maxiter = 20, double abs_diff = 1e-6, int verbose = 0){

	int n = n_x * n_y * n_z;	double old_val = 0.0;	int bad_count_n_o = 0, bad_count_o = 0;
	Vector_eig r2 = Vector_eig::Ones(3);
	EM_opt<float> f(r2);
	Vector_eig lb(3*n+6+1+1), ub(3*n+6+1+1); 
	lb = Vector_eig::Zero(3*n+6+1+1, 1); 	ub = Vector_eig::Ones(3*n+6+1+1, 1);	ub.segment(0, n) *= 255;
	
	
	
	
	// Bound on Cholesky etc:
	ub.segment(3*n, 6+1+1) *= 255;
	lb.segment(3*n, 6+1+1) = -255*Vector_eig::Ones(6+1+1);
	lb(3*n) = 1e-5; lb(3*n+3) = 1e-5;lb(3*n+5) = 1e-5;	lb(3*n+6) = 1e-5;lb(3*n+7) = 1e-5;
	
	//ub.segment(3*n, 6+1+1) *= 15;		// Get's sqaured
	//lb.segment(3*n, 6+1+1) = 0.00001*Vector_eig::Ones(6);
	
	
	
	f.r = r;	f.TE = TE_example;		f.TR = TR_example;
	f.setLowerBound(lb);	f.setUpperBound(ub);
	f.n_x = n_x; f.n_y = n_y; f.n_z = n_z;
	f.sigma = sigma;
	
	// For the extra check:
	f.ub = ub; 		f.lb = lb;
	
	double beta_x = 0.1, beta_y = 0.1, beta_z = 0.1;			// beta_z -- ??
	Psi_inv = Eigen::Matrix3f::Identity()*1.1;	// 1.1 is just for test of ub
	
	Matrix_eig W_old = W_init;
	f.W_old = W_old;
	Vector_eig param_new = to_param_vec(W_init, Psi_inv, beta_x, beta_y);		//W_init.row(i);
	Vector_eig param_old = param_new;
	Debug0 ("f(x) at first: \n");
	Debug0(f.value(param_new));
	
	int iter = 0;
	double best_val;
	
	// EM loop - M step//
	while(iter < maxiter){
	
		Debug1("\n--------------------------------------------------------------------------\nIteration: " << iter++ << "\n");
		
		W_old.col(0) = param_old.segment(0,  n);
		W_old.col(1) = param_old.segment(n,  n);
		W_old.col(2) = param_old.segment(2*n,n);
		f.W_old = W_old;
		old_val = f.value(param_old);
		
		// Track the best:
		float current_best_val = 1.0e+15;
		
		
		//Solve:
		cppoptlib::LbfgsbSolver<EM_opt<float>> solver;
		solver.minimize(f, param_new);
		Debug1("Solved for this iteration!\n");
		
		// Track the best:
		param_new = f.current_best_param;
		best_val = f.current_best_val;
		
		std::cout << "f(param_new) in argmin: " << best_val << "\t old val:" << old_val << "\n";
		std::cout << "Solver status: " << solver.status() << "\n";	//Guess: bad reports: under constraints => grad is not ~0 
		std::cout << "Final criteria values: " << "\n" << solver.criteria() << "\n";
		
		
		//moving on part:
		if(abs_sum(param_old-param_new) <= abs_diff){
			std::cout << "Stopped after " << iter << " iterations" << "\n";
			break;
		}
		
		if(best_val>=old_val){ 
			Debug1("Value not decreased!! old val: " << old_val << ";\t new val: " << best_val);
			bad_count_o++;
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
	Matrix3f_eig L_mat = to_L_mat(temp_L);
	Psi_inv = from_Cholesky(L_mat);
	
	//beta
	beta(0) = param_new(3*n+6); beta(1) = param_new(3*n+7); beta(2) = beta_z;
}


// abs -- done



/***************************************************
**************** Information Matrix ****************
****************************************************/


SpMat Hessian_mat_without_MRF(Matrix_eig W, Matrix3f_eig Psi_inv, Vector_eig beta, 
                   Vector_eig TE, Vector_eig TR, Vector_eig sigma, Matrix_eig r, 
                   int n_x, int n_y, int n_z){
	
	auto time_1_hess = std::chrono::high_resolution_clock::now();
	Debug1("Hessian calculation started without MRF");
	Matrix_eig v = v_mat(W, TE, TR);
	int n = n_x * n_y * n_z;	//n
	int m = v.cols();			//m
	double temp = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;
	SpMat Gamma_inv = Lambda(beta, n_x, n_y, n_z);
	SpMat W_hess(3*n, 3*n);
	W_hess.reserve( VectorXi::Constant(3*n, 9) );		// Reserve 9 non-zero's per column
	Debug1("Hessian calculation allocated");
	
	// W - grad part//
	int i = 0, k = 0, k1 = 0, j = 0;
	Vector_eig temp_vec(3);
	for(i = 0; i < n; ++i) {
		temp_vec = W.row(i);
		
		for(k = 0; k < 3; ++k) {
			for(k1 = 0; k1 < 3; ++k1) {
				temp = 0.;									// Missed this -- correct this.
				for(j = 0; j < m ; ++j) {
					tmp2 = r(i,j)/SQ(sigma(j));
					tmp3 = -v(i,j)/SQ(sigma(j)) + tmp2*besselI1_I0(tmp2*v(i,j));
					temp += tmp3* simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(temp_vec, TE, TR, j, k, k1);
					
					tmp2 *= v(i,j);
					tmp3 = (1 + ratio_bessel_20(tmp2) - 2*SQ(besselI1_I0(tmp2)) );
					tmp2 /= v(i,j);
					tmp4 = -1/SQ(sigma(j)) + 0.5*SQ(tmp2)*tmp3;
					temp += tmp4 * simple_dee_v_ij_dee_W_ik(temp_vec, TE, TR, j, k) * simple_dee_v_ij_dee_W_ik(temp_vec, TE, TR, j, k1);
				}
				W_hess.insert(i+k*n, i+k1*n) = temp;
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







int main(int argc, char * argv[]) {
	
	if (argc != 4) {
		fprintf(stderr, "\nUsage: %s <file_name> <SD_file_name> <will_write_to_a_file?> <temp_val> \n", argv[0]);
		exit(EXIT_FAILURE);
	}
	char *data_file, *sd_file;
	data_file = argv[1]; 	sd_file = argv[2]; 	char will_write = *(argv[3])-48;		// Converted from ascii
	short our_dim[8];

	Vector_eig TE_example((Vector_eig(12) << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1).finished());
	Vector_eig TR_example((Vector_eig(12) << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3).finished());


	// Test:
	// Main step:
	Matrix_eig r = Preprocess_data(data_file, our_dim, will_write);
	Debug2("Preprocessing done");
	
	int do_least_sq = 0;	// 0 Subrata -- least sq have better initial likelihood-but stucks and gives nan in some value
	if(do_least_sq){
		Debug0("Doing least sq estimate!");
	}
	Matrix_eig W_init = Init_val(r, TE_example, TR_example, our_dim, do_least_sq, will_write);
	Debug2("W initial done");
	check_nan(W_init);
	show_head(W_init);
	
	Vector_eig sigma = read_sd(sd_file, our_dim[4]);
	Debug0("sigma: " << sigma.transpose());

	int n = our_dim[1]*our_dim[2]*our_dim[3];
	Eigen::Matrix3f Psi_inv_init = Eigen::Matrix3f::Identity();
	Vector_eig beta_init = Vector_eig::Zero(3);
	EM_solve(W_init, Psi_inv_init, beta_init, TE_example, TR_example, r, our_dim[1], our_dim[2], our_dim[3], sigma);
	show_head(W_init);



	Hessian_mat_without_MRF(W_init, Psi_inv_init, beta_init, TE_example,TR_example, sigma, r, our_dim[1], our_dim[2], our_dim[3]);

	return 0;
}


