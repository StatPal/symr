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


/**
*
* File to get initial values through Least Squares :
* Example use in main :


int main(int argc, char * argv[]){
	if (argc != 4) {
		fprintf(stderr, "\nUsage: %s <file_name> <SD_file_name> <will_write_to_a_file?> <temp_val> \n", argv[0]);
		exit(EXIT_FAILURE);
	}
	char *data_file, *sd_file;
	data_file = argv[1];	sd_file = argv[2];	char will_write = *(argv[3])-48;		// Converted from ascii
	short our_dim[8];
	// Values of TE and TR (in seconds)
	Vector_eig TE_example((Vector_eig(12) << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1).finished());
	Vector_eig TR_example((Vector_eig(12) << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3).finished());
	
	// Main step:
	int do_least_sq = 1;
	Matrix_eig_row r = Preprocess_data(data_file, our_dim, will_write);
	Matrix_eig_row W_init = Init_val(r, TE_example, TR_example, our_dim, do_least_sq, will_write);
	
	return 0;
}

*/



#ifndef INIT_VAL
#define INIT_VAL



#include "functions_gen.hpp"
#include "read_files.hpp"



// #include "../CppNumericalSolvers/include/cppoptlib/meta.h"
#include "../CppNumericalSolvers/include/cppoptlib/boundedproblem.h"
#include "../CppNumericalSolvers/include/cppoptlib/solver/lbfgsbsolver.h"









/**
* Optimization template for Least Square Estimates of the parameters.
*/
template<typename T>
class Least_Sq_est : public cppoptlib::BoundedProblem<T> {
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	typedef Matrix_eig_row TMatrix_row;
	
	TVector r_row;
	double fx;


  public:
	Least_Sq_est() : 
		cppoptlib::BoundedProblem<T>(3){}
	


	TVector TE, TR, lb, ub, v_new;
	int i;
	
	
	void update_size(){
		v_new = TVector::Zero(TE.size());
	}
	
	

	// Track the best:
	Vector_eig current_best_param;
	double current_best_val = 1.0e+15;


	// Objective function, to be minimized:
	T value(const TVector &x) {
		Bloch_vec(x, TE, TR, v_new);
		fx = (r_row - v_new).squaredNorm();
		
		// Track the best:
		if(fx < current_best_val){
			if(check_bounds_vec_3(x, lb, ub) == 0){
				current_best_param = x;
				current_best_val = fx;
			}
		}
				
		return (fx);
	}


	// grad of the value: 
	void gradient(const TVector &x, TVector &grad) {
		Bloch_vec(x, TE, TR, v_new);
		grad << 0,0,0;
		int m = TR.size();
		
		
		for(int j = 0; j < m; ++j){
		
			grad[0] -= 2*(r_row(j) - v_new(j)) * 
						std::exp(TE(j)*std::log(x(2))) * 
						(1-std::exp(TR(j)*std::log(x(1))));
			
			grad[1] -= 2*(v_new(j) - r_row(j)) * 
						x(0)*TR(j) * std::exp(TE(j)*log(x(2))) * 
						std::exp((TR(j)-1)*std::log(x(1)));
			
			grad[2] -= 2*(r_row(j) - v_new(j)) * 
						x(0)*TE(j) * std::exp((TE(j)-1)*log(x(2))) * 
						(1-std::exp(TR(j)*std::log(x(1))));
		}
		if(i == 2)
			DebugLS("grad: " << grad.transpose() << "\n" );
	}


};




/**
* Least square solution
* 	Input:  
*			W, TE_example, TR_example, r, r_scale, TE_scale, TR_scale
* 			
* 			W is changed
*/
void least_sq_solve(Matrix_eig_row &W, 
					const Vector_eig &TE_example, const Vector_eig &TR_example, 
                    const Matrix_eig_row &r, double r_scale, double TE_scale, double TR_scale){


	Debug0("Doing Least Square Estimate!");
	auto time_1_lsq = std::chrono::high_resolution_clock::now();
	
	Least_Sq_est<double> f;
	Vector_eig x(3), lb(3), ub(3); 
	
	//Bounds of rho, W1, W2:
	lb << 0.0001, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	for(int i = 1; i < 3; ++i){
		if(lb[i]<1.0e-8){
			lb[i] = 1.0e-8;
		}
	}
	Debug1("lb inside LS: " << lb.transpose());
	Debug1("ub inside LS: " << ub.transpose() << "\n");
	
	
	f.TE.noalias() = TE_example;	f.TR.noalias() = TR_example;
	f.setLowerBound(lb);	f.setUpperBound(ub);		f.lb.noalias() = lb; 	f.ub.noalias() = ub;
	f.update_size();
	
	double old_val = 0.0, fx;
	int n = r.rows(), bad_count_o = 0, bad_count_o_2 = 0, bad_bound_1 = 0, bad_bound_2 = 0, nan_count = 0;



	
	
	// Declaring the solver
	cppoptlib::LbfgsbSolver<Least_Sq_est<double>> solver;
	cppoptlib::Criteria<double> crit_LS = cppoptlib::Criteria<double>::defaults();
	crit_LS.iterations = 100;
	solver.setStopCriteria(crit_LS);
	
	
	
	
	// See https://bisqwit.iki.fi/story/howto/openmp/#PrivateFirstprivateAndSharedClauses for modifications also
	
	// Loop of 
	#pragma omp parallel for default(none) firstprivate(f, solver) private (x, old_val, fx)  shared(W, bad_count_o, nan_count, bad_count_o_2, r, TE_example, TR_example, n, Rcpp::Rcout)
	for(int i = 0; i < n; ++i){
	
		if(i % 100000 == 0){
			Debug1("i: "<< i);
		}
		
		// Track the best:
		f.current_best_val = 1.0e+15;
		
		f.i = i;
		x = W.row(i);
		f.r_row = r.row(i);
		
	
		//Print initial values:
		DebugLS ("value of i: " << i << "\t x at first: " << x.transpose());
		old_val = f.value(x);
		DebugLS ("f(x) at first: \n" << old_val) ;
		
	
		//Solve:
		solver.minimize(f, x);
		fx = f.value(x);
		
		// Track the best:
		x = f.current_best_param;
		fx = f.current_best_val;
		DebugLS("f(param_new) in argmin: " << fx << "\t x:" << x.transpose());
		
		
		
		if(fx>=old_val){
			Debug2("Value not decreased!! old x:" << W.row(i) << " val: " << old_val << 
						";\t x: " << x.transpose() << " val: " << fx << " i:" << i);
			bad_count_o++;
			if(fx>old_val){
				bad_count_o_2++;
			}
		} else {
			if(check_nan_vec(x) == 0){				// Added later, to catch NaN - Subrata
				W.row(i) = x;
			} else {
				Debug0("nan in LS estimate. \n" << "i: " << i << ", x: " << x.transpose());
				nan_count++;
			}
		}
		
	}
	
	Debug0("Number of bad cases in Initial value determination:" << bad_count_o << 
			" and worse: " << bad_count_o_2 << 
			" and bad init bounds:" << bad_bound_1 << " and " << bad_bound_2);
	
	auto time_2_lsq = std::chrono::high_resolution_clock::now();
	auto duration_lsq = std::chrono::duration_cast<std::chrono::microseconds>(time_2_lsq - time_1_lsq);
	Debug1("Time taken total Least Square: " << duration_lsq.count() << " microseconds\n");
}






/**
* Read the data and
* replace the zeros by small number (0.5)
*/
Matrix_eig_row Preprocess_data(char* const data_file, short our_dim[8], char will_write = 0){
	
	// Load the data file first. 
	Matrix_eig_row r = Read_nift1(data_file, our_dim, will_write);
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





/**
* Creates the initial matrix W
* if do_least_sq is 1, 
* 	it gives the least square solution.
*/
Matrix_eig_row Init_val(const Matrix_eig_row &r, 
                        const Vector_eig &TE_example, const Vector_eig &TR_example, 
                        short our_dim[8], 
                        double r_scale, double TE_scale, double TR_scale, 
                        double W_1_init = exp(-1/2.0), double W_2_init = exp(-1/0.1), 
                        int do_least_sq = 1, char will_write = 0){


	//Primary Initial value for test//
	int n = our_dim[1]*our_dim[2]*our_dim[3];
	Matrix_eig_row W = Matrix_eig_row::Ones(n, 3);
	
	
	// Ad hoc initial values:
	W.col(0) = r.rowwise().mean().transpose();
	W.col(1) *= W_1_init;
	W.col(2) *= W_2_init;
	for(int i = 0; i < r.rows(); ++i){
		if(W(i, 0)>450.0){
			W(i, 0) = 425.0;
		}
		if(W(i, 0) < 0.0001){		// Added later to recover from sub lb[0] case possiblw due to scaling down
			W(i, 0) = 0.0001;
		}

		
	}
	Debug1("1st level preprocess of initial value done!\n----------------\n----------------\n");



	if(DEBUG_ANOTHER_LEVEL){
		show_head(W);
	}
	if(do_least_sq){
		least_sq_solve(W, TE_example, TR_example, r, r_scale, TE_scale, TR_scale);
	}
	if(DEBUG_ANOTHER_LEVEL){
		std::cout << "After the operation:";
		show_head(W);
	}
	
	Debug1("Initial value done!\n----------------\n");

	return W;
}




/**
* Performance matrices:
* W: 		parameter matrix:
* test: 	test set image matrix
* TE_test:			
* TR_test:			
* sigma_test:		
* v_type 
	1 -- compared with v
	2 -- compared with mode of rice distn (not implemented yet - need another set of optimization)
	3 -- compared with rice mean
	
* measure_type: 
  	1 -- abs deviation from mean(not median?)
  	2 -- squared deviation from mean
  
* Scale: Scaled measure vs Unscaled measure 
  	0 -- no
  	1 -- yes
		 Scaling factor is sd / MAD w.r.t. mean
*/
Vector_eig Performance_test(const Matrix_eig_row &W, const Matrix_eig_row &test, 
							const Vector_eig &TE_test, const Vector_eig &TR_test, const Vector_eig &sigma_test, 							const Eigen::Matrix<char, Eigen::Dynamic, 1> &black_list,
							int v_type = 1, int measure_type = 1, int scale = 1, int verbose = 0){


	
	int n_test = TE_test.size(), n = W.rows();
	assert(sigma_test.size() == n_test);
	Vector_eig Performance_test = Vector_eig::Zero(n_test);
	Matrix_eig_row Perf_mat = Matrix_eig_row::Zero(W.rows(), n_test);		// just testing 
	
	Vector_eig tmp(n_test);
	
	Vector_eig v_new = Vector_eig::Zero(n_test);
	Vector_eig v_star(n_test);
	
	int fg_num = 0;
	
	
	// Not exactly correct: Subrata - Check
	
	// #pragma omp parallel for default(none) firstprivate(v_new, v_star, tmp) shared(W, n, test, n_test, TE_test, TR_test, sigma_test, v_type, measure_type, verbose, std::cout, Perf_mat)		// reduction(+:Performance_test)
	for(int i = 0; i < n; ++i) {
		if(black_list(i) == 0){
			fg_num++;
			Bloch_vec(W.row(i), TE_test, TR_test, v_new);			// v_{ij}
			
			
			if(v_type == 1){
				v_star = v_new;
			} else if (v_type == 2){
				// Need to be done - mode of rice distn to be calculated with NR method
			} else if (v_type == 3){
				for(int j = 0; j < n_test; ++j){
					v_star(j) = mean_rice(v_new(j), sigma_test(j));
				}
			}
			
			tmp = (v_star - test.row(i).transpose()).array().abs();
		
			if(measure_type == 2){
				for(int j = 0; j < n_test; ++j){
					tmp(j) = SQ(tmp(j));
				}
			}
			
			
			if(verbose){
				if(i < 100){
					Debug1("i: " << i << ", v_new:" << v_new.transpose() <<  ", v_star:" << v_star.transpose() << 
							",\n test: " << test.row(i) <<  "\n v_star - test_row" << v_star.transpose() - test.row(i) << 
							 ", tmp: " << tmp.transpose() << "\n");
				}			
			}
			
			int bad_ct = 0;
			for(int j = 0; j < tmp.size(); ++j){
				if(std::isnan(tmp(j))){
					bad_ct++;
					Debug1("i: " << i << ", W.row(i): " << W.row(i) << 
							"\n, v_new:" << v_new.transpose() <<  ", v_star:" << v_star.transpose() << 
							",\n test: " << test.row(i) <<  "\n v_star - test_row" << v_star.transpose() - test.row(i) << 
							 ", tmp: " << tmp.transpose() << "\n");
					//break;
					exit(EXIT_FAILURE);
				}
			}
			
			
			Perf_mat.row(i) = tmp;
			// Performance_test = Performance_test + tmp;				// This is main
		}
	}
	
	// Performance_test = Performance_test/W.rows();
	Performance_test = Perf_mat.array().colwise().sum()/fg_num;
	
	
	
	if(verbose){
		 Debug0("Performance_test: " << Performance_test.transpose());
	}
	
	if(measure_type == 2){
		for(int j = 0; j < n_test; ++j){
			Performance_test[j] = std::sqrt(Performance_test[j]);
		}
	}
	
	if(verbose){
		 Debug0("Performance_test: " << Performance_test.transpose());
	}
	
	
	if(scale){
		double scale_factor = 1.0;
		for(int j = 0; j < n_test; ++j){
			if(measure_type == 1){
				scale_factor = abs_dev_mean(test.col(j)); 
			} else if (measure_type == 2){
				scale_factor = std::sqrt(var(test.col(j)));
			}
			Performance_test[j] /= scale_factor;
		}
	}
	
	
	return Performance_test;
}




#endif	/* INIT_VAL */
