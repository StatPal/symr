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
	Eigen::VectorXf TE_example((Eigen::VectorXf(12) << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1).finished());
	Eigen::VectorXf TR_example((Eigen::VectorXf(12) << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3).finished());
	
	// Main step:
	int do_least_sq = 1;
	Eigen::MatrixXf r = Preprocess_data(data_file, our_dim, will_write);
	MatrixXf W_init = Init_val(r, TE_example, TR_example, our_dim, do_least_sq, will_write);
	
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
	Eigen::VectorXf current_best_param;
	float current_best_val = 1.0e+15;


	T value(const TVector &x) {
		TVector v_new = Bloch_vec(x, TE, TR);
		Debug2("x: " << x.transpose() << "; value: " << (r.row(i).transpose() - v_new).squaredNorm()); 
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
		grad[0] = 0; grad[1] = 0; grad[2] = 0;
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
		Debug2("grad: " << grad.transpose() << "\n" );
	}
};


void least_sq_solve(Eigen::MatrixXf &W, Eigen::VectorXf TE_example, Eigen::VectorXf TR_example, Eigen::MatrixXf &r){

	auto time_1_lsq = std::chrono::high_resolution_clock::now();
	
	Eigen::VectorXf r2 = Eigen::VectorXf::Ones(3) * 0.5;
	Least_Sq_est<float> f(r2);
	Eigen::VectorXf x(3), lb(3), ub(3); 
	lb << 0, 0, 0;
	ub << 255, 1, 1;
	f.r = r;	f.TE = TE_example;	f.TR = TR_example;	f.setLowerBound(lb);	f.setUpperBound(ub);
	double old_val = 0.0;
	int bad_count_o = 0;
	// Declaring the solver
	cppoptlib::LbfgsbSolver<Least_Sq_est<float>> solver;

	for(int i = 0; i < r.rows()/10000; ++i){
		
		if(i==100000 || i==200000 || i==300000 || i==400000 || i==500000 || i==600000 || i==700000 || i==800000 || i==900000 ){
			Debug1("i: "<< i);
		}
		
		
		x = W.row(i);
		f.i = i;
		// Track the best:
		float current_best_val = 1.0e+15;
	
		//Print initial values:
		Debug2 ("value of i: " << i << "\t x at first: " << x.transpose());
		Debug2 ("f(x) at first: \n" << f.value(x)) ;
		old_val = f.value(x);
	
		//Solve:
		solver.minimize(f, x);
		Debug2("argmin: " << x.transpose() << ";\tf(x) in argmin: " << f(x)) ;
		Debug2("Solver status: " << solver.status() );	//Guess: bad reports: under constraints => grad is not ~0 
		Debug2("Final criteria values: " << "\n" << solver.criteria());
	
		// Track the best:
		x = f.current_best_param;
		double fx = f.current_best_val;
		Debug2("f(param_new) in argmin: " << fx << "\t old val:" << old_val);
		Debug2(niter << " iterations");
		
		
		if(fx>=old_val){
			Debug2("Value not decreased!! old x:" << W.row(i) << " val: " << old_val << ";\t x: " << x.transpose() << " val: " << fx << " i:" << i);
			bad_count_o++;
		} else {
			W.row(i) = x;
		}
		
		
		f.current_best_val = 1.0e+15;
	}
	
	Debug1("Number of bad cases in Initial value determination:" << bad_count_o);
	
	auto time_2_lsq = std::chrono::high_resolution_clock::now();
	auto duration_lsq = std::chrono::duration_cast<std::chrono::seconds>(time_2_lsq - time_1_lsq);
	Debug1("Time taken total loop: " << duration_lsq.count() << " seconds\n");
}


Eigen::MatrixXf Preprocess_data(char* const data_file, short our_dim[8], char will_write = 0){
	
	// Load the data file first. 
	Eigen::MatrixXf r = Read_nift1(data_file, our_dim, will_write);
	for(int i = 0; i < our_dim[1]*our_dim[2]*our_dim[3]; ++i){
		for(int j = 0; j < our_dim[4]; ++j){
			if(r(i, j) == 0){
				r(i, j) = 0.5;			// Just added a small value to remove -Inf problem of likelihood.
			}
		}
	}
	return r;
}


Eigen::MatrixXf Init_val(Eigen::MatrixXf r, Eigen::VectorXf TE_example, Eigen::VectorXf TR_example, short our_dim[8], int do_least_sq = 1, char will_write = 0){

	//Primary Initial value for test//
	int n = our_dim[1]*our_dim[2]*our_dim[3];
	Matrix_eig W = Matrix_eig::Ones(n, 3);
	W.col(0) = r.rowwise().mean().transpose();
	W.col(2) *= 0.9;
	for(int i = 0; i < r.rows(); ++i){
		if(W(i, 0) == 0.5){
			W(i, 1) = 0.2;
		} else {
			W(i, 0) *= 1.2;
			W(i, 1) = 0.8;
		}
	}
	
	// See Init_value_3 and Init_value_4 for more explanations.
	// If same -- 1-1.1, 0.1 - 0.3, 0.9
	// Not same -- 1, 0.7-0.8, 0.9
	Debug1("1st level preprocess of initial value done!\n----------------\n");


	if(DEBUG_ANOTHER_LEVEL){
		show_head(W);
	}
	if(do_least_sq){
		least_sq_solve(W, TE_example, TR_example, r);
	}
	if(DEBUG_ANOTHER_LEVEL){
		std::cout << "After the operation:";
		show_head(W);
	}
	
	Debug1("Initial value done!\n----------------\n");

	return W;
}




#endif	/* INIT_VAL */
