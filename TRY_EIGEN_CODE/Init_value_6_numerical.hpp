/**
*
* File to get initial values through least squares:
* Example use in main:


int main(int argc, char * argv[]){
	if (argc != 4) {
		fprintf(stderr, "\nUsage: %s <file_name> <SD_file_name> <will_write_to_a_file?> <temp_val> \n", argv[0]);
		exit(EXIT_FAILURE);
	}
	char *data_file, *sd_file;
	data_file = argv[1];	sd_file = argv[2];	char will_write = *(argv[3])-48;		// Converted from ascii
	short our_dim[8];
	// Values of TE and TR (in seconds)
	Eigen::VectorXd TE_example((Eigen::VectorXd(12) << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1).finished());
	Eigen::VectorXd TR_example((Eigen::VectorXd(12) << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3).finished());
	
	// Main step:
	int do_least_sq = 1;
	Eigen::MatrixXd r = Preprocess_data(data_file, our_dim, will_write);
	MatrixXd W_init = Init_val(r, TE_example, TR_example, our_dim, do_least_sq, will_write);
	
	return 0;
}

*/





#include "scheme_new_numerical.hpp"
#include "Read_files_2.hpp"


#include "../optim_cpp_solver/include/cppoptlib/meta.h"
#include "../optim_cpp_solver/include/cppoptlib/boundedproblem.h"
#include "../optim_cpp_solver/include/cppoptlib/solver/lbfgsbsolver.h"


#ifndef INIT_VAL
#define INIT_VAL

//const Vector_eig TE_example((Vector_eig(12) << 10, 15, 20, 10, 30, 40, 10, 40, 80, 10, 60, 100).finished());
//const Vector_eig TR_example((Vector_eig(12) << 600, 600, 600, 1000, 1000, 1000, 2000, 2000, 2000, 3000, 3000, 3000).finished());

//const Vector_eig TE_example((Vector_eig(12) << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1).finished());
//const Vector_eig TR_example((Vector_eig(12) << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3).finished());





template<typename T>
class Least_Sq_est : public cppoptlib::BoundedProblem<T> {
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	TMatrix r;
	TVector r2;

  public:
	Least_Sq_est(const TVector y_) :											// Is it constructor?
		cppoptlib::BoundedProblem<T>(y_.size()), r2(y_){}
	// What happens to r(y_) when r is a matrix. 
	// r(y_){} changed to r2(y_){}													// Know more about it.

	TVector TE, TR;
	int i;
	
	// Track the best:
	Eigen::VectorXd current_best_param;
	double current_best_val = 1.0e+15;


	T value(const TVector &x) {
		TVector v_new = Bloch_vec(x, TE, TR);
		Debug3("x: " << x.transpose() << "; value: " << (r.row(i).transpose() - v_new).squaredNorm()); 
		double fx = (r.row(i).transpose() - v_new).squaredNorm();
		
		// Track the best:
		if(fx < current_best_val){
			current_best_param = x;
			current_best_val = fx;
		}
		
		return (fx);
	}


	void gradient(const TVector &x, TVector &grad) {
		TVector v_new = Bloch_vec(x, TE, TR);
		grad << 0,0,0;						//grad[0] = 0; grad[1] = 0; grad[2] = 0;
		int m = TR.size();
		
		if(x(1)==0){														// To avoid numerical problem at boundary.
			for(int j = 0; j < m; ++j){
				grad[0] -= 2*(r(i, j) - v_new(j)) * exp(TE(j)*log(x(2))) ;
				grad[1] -= 0;
				grad[2] -= 2*(r(i, j) - v_new(j)) * x(0)*TE(j) * exp((TE(j)-1)*log(x(2))) ;
			}
		} else { 
			for(int j = 0; j < m; ++j){
				grad[0] -= 2*(r(i, j) - v_new(j)) * exp(TE(j)*log(x(2))) * (1-exp(TR(j)*log(x(1))));
				grad[1] -= 2*(v_new(j) - r(i, j)) * x(0)*TR(j) * exp(TE(j)*log(x(2))) * exp((TR(j)-1)*log(x(1)));
				grad[2] -= 2*(r(i, j) - v_new(j)) * x(0)*TE(j) * exp((TE(j)-1)*log(x(2))) * (1-exp(TR(j)*log(x(1))));
			}
		}
		Debug3("grad: " << grad.transpose() << "\n" );
	}


};


void least_sq_solve(Eigen::MatrixXd &W, Eigen::VectorXd TE_example, Eigen::VectorXd TR_example, Eigen::MatrixXd &r, 
                    double TE_scale, double TR_scale){


	Debug0("Doing Least Square Estimate!");
	auto time_1_lsq = std::chrono::high_resolution_clock::now();
	
	Eigen::VectorXd r2 = Eigen::VectorXd::Ones(3) * 0.5;
	Least_Sq_est<double> f(r2);
	Eigen::VectorXd x(3), lb(3), ub(3); 
	
	//Bounds of rho, W1, W2:
	lb << 0.0, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	
	
	
	f.r = r;	f.TE = TE_example;	f.TR = TR_example;
	f.setLowerBound(lb);	f.setUpperBound(ub);
	double old_val = 0.0;
	int bad_count_o = 0, bad_count_o_2 = 0, bad_bound_1 = 0, bad_bound_2 = 0;



	// Declaring the solver
	cppoptlib::LbfgsbSolver<Least_Sq_est<double>> solver;
	
	for(int i = 0; i < r.rows()/1; ++i){
	
		if(i==100000 || i==200000 || i==300000 || i==400000 || i==500000 || i==600000 || i==700000 || i==800000 || i==900000 ){
			Debug1("i: "<< i);
		}
		
		
		x = W.row(i);
		
		//Check Bounds:
		for(int j = 0; j < 3; ++j){
			if(x[0]<lb[0] || x[1]<lb[1] || x[2]<lb[2]){
				Debug1("Crossed lower bound initially!");
				bad_bound_1++;
			}
			if(x[0]>ub[0] || x[1]>ub[1] || x[2]>ub[2]){
				Debug1("Crossed upper Bound initially!");
				Debug1("x " << x.transpose() << " ub: " << ub.transpose());
				bad_bound_2++;
			}
		}
		
		
		f.i = i;
		// Track the best:
		double current_best_val = 1.0e+15;
	
		//Print initial values:
		Debug3 ("value of i: " << i << "\t x at first: " << x.transpose());
		Debug3 ("f(x) at first: \n" << f.value(x)) ;
		old_val = f.value(x);
	
		//Solve:
		solver.minimize(f, x);
		Debug3("argmin: " << x.transpose() << ";\tf(x) in argmin: " << f(x)) ;
		Debug3("Solver status: " << solver.status() );	//Guess: bad reports: under constraints => grad is not ~0 
		Debug3("Final criteria values: " << "\n" << solver.criteria());
	
		// Track the best:
		x = f.current_best_param;
		double fx = f.current_best_val;
		Debug3("f(param_new) in argmin: " << fx << "\t old val:" << old_val);
		//Debug2(niter << " iterations");			// Wait! niter not defined ????
		
		
		if(fx>=old_val){
			Debug2("Value not decreased!! old x:" << W.row(i) << " val: " << old_val << ";\t x: " << x.transpose() << " val: " << fx << " i:" << i);
			bad_count_o++;
			if(fx>old_val){
				bad_count_o_2++;
			}
		} else {
			W.row(i) = x;
		}
		
		
		f.current_best_val = 1.0e+15;
	}
	
	Debug0("Number of bad cases in Initial value determination:" << bad_count_o << " and worse: " << bad_count_o_2 << " and bad init bounds:" << bad_bound_1 << " and " << bad_bound_2);
	
	auto time_2_lsq = std::chrono::high_resolution_clock::now();
	auto duration_lsq = std::chrono::duration_cast<std::chrono::seconds>(time_2_lsq - time_1_lsq);
	Debug1("Time taken total Least Square: " << duration_lsq.count() << " seconds\n");
}


Eigen::MatrixXd Preprocess_data(char* const data_file, short our_dim[8], char will_write = 0){
	
	// Load the data file first. 
	Eigen::MatrixXd r = Read_nift1(data_file, our_dim, will_write);
	for(int i = 0; i < our_dim[1]*our_dim[2]*our_dim[3]; ++i){
		for(int j = 0; j < our_dim[4]; ++j){
			if(r(i, j) == 0){
				r(i, j) = 0.5;			// Just added a small value to remove -Inf problem of likelihood.
			}
		}
	}
	Debug1("Preprocessing done!");
	return r;
}



Eigen::MatrixXd Init_val(Eigen::MatrixXd r, Eigen::VectorXd TE_example, Eigen::VectorXd TR_example, short our_dim[8], 
                         double TE_scale, double TR_scale, double W_1_init = exp(-1/2.0), double W_2_init = exp(-1/0.1), 
                         int do_least_sq = 1, char will_write = 0){

	//Primary Initial value for test//
	int n = our_dim[1]*our_dim[2]*our_dim[3];
	Matrix_eig W = Matrix_eig::Ones(n, 3);
	//show_dim(W);
	//show_dim(r);
	
	// Ad hoc initial values:
	W.col(0) = r.rowwise().mean().transpose();
	W.col(1) *= W_1_init;
	W.col(2) *= W_2_init;
	for(int i = 0; i < r.rows(); ++i){
		if(W(i, 0)>450){
			W(i, 0) = 425;
		}
	}
	Debug1("1st level preprocess of initial value done!\n----------------\n----------------\n");



	if(DEBUG_ANOTHER_LEVEL){
		show_head(W);
	}
	if(do_least_sq){
		least_sq_solve(W, TE_example, TR_example, r, TE_scale, TR_scale);
	}
	if(DEBUG_ANOTHER_LEVEL){
		std::cout << "After the operation:";
		show_head(W);
	}
	
	Debug1("Initial value done!\n----------------\n");

	return W;
}



Vector_eig Performance_test(Matrix_eig W, Matrix_eig test, Vector_eig TE_test, Vector_eig TR_test,
							Vector_eig sigma_test, 
							int v_type = 1, int measure_type = 1){

	int n_test = TE_test.size();
	assert(sigma_test.size() == n_test);
	Vector_eig Performance_test = Vector_eig::Zero(n_test);
	Matrix_eig Perf_mat = Matrix_eig::Zero(W.rows(), n_test);		// test 
	for(int i = 0; i < W.rows(); ++i){
		Vector_eig v_new = Bloch_vec(W.row(i), TE_test, TR_test);			// v_{ij}
		Vector_eig v_star(n_test);
		if(v_type == 1){
			v_star = v_new;
		} else if (v_type == 2){
			// Need to be done - mode of rice distn to be calculated with NR method
		} else if (v_type == 3){
			for(int j = 1; j < n_test; j++){
				v_star(j) = mean_rice(v_new(j), sigma_test(j));
			}
		}
		
		Vector_eig tmp = (v_star - test.row(i).transpose()).array().abs();
		
		if(measure_type == 2){
			for(int j = 0; j < n_test; ++j){
				tmp(j) = SQ(tmp(j));
			}
		}
		
		
		//std::cout << tmp.transpose() << "\n";
		Perf_mat.row(i) = v_star.transpose() - test.row(i);
		
		
		Performance_test = Performance_test + tmp;
	}
	std::cout << Perf_mat.colwise().mean() << " and " <<  Perf_mat.array().abs().colwise().mean() << "\n";
	//Performance_test = Perf_mat.array().abs().colwise().mean();
	
	
	Performance_test = Performance_test/W.rows();
	
	if(measure_type == 2){
		for(int j = 0; j < n_test; ++j){
			Performance_test[j] = std::sqrt(Performance_test[j]);
		}
	}
	
	return Performance_test;
}
// Do scaled version:
// For scaled version, RMSPE is rescaled using crude SD of j-th image
// MAPE is rescaled using crude MAD about the mean(??? isn't it better to take median) of j-th image



#endif	/* INIT_VAL */
