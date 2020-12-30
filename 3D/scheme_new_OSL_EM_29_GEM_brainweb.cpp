/**
* 
* Multicycle EM (AECM) algorithm
* Psi and beta are updated at every loop
* Checkerboard structure implemented.
* Some more advancement in MRF likeli num for faster AECM
* E step is done seperately in that class of optimization. 

* Parallelized


* To compile:

g++ scheme_new_OSL_EM_29_GEM_brainweb.cpp -o test_29_3D_brainweb -I ~/program/eigen3 -O3 -lgsl -lgslcblas -lm -fopenmp -DEIGEN_DONT_PARALLELIZE




./test_29_3D_brainweb ../Read_Data/brainweb_all.nii Dummy_sd_brainweb.txt 0

./test_29_3D_brainweb ../data/brainweb_all.nii Dummy_sd_brainweb.txt 0

nohup ./test_29_3D_brainweb ../Read_Data/brainweb_all.nii Dummy_sd_brainweb.txt 0 > test_29_3D_brainweb.out & 




Changes:

Black listed pixels


MRF estimation in a new way


Rewriting everything in a new way for parallel


TODO:
Write the update neighbour step inside the Likeli_optim - or something like that

Write E step before M step - would possibly have easier to compare in difference of W's
or
Have another W_old just to compare: W_old_reserve 
----------WAIT, but in this case, W_old and W are always same - what's the utility?????

* 
*/




#include "functions_gen.hpp"
#include "read_files.hpp"
#include "functions_LS_and_init_value.hpp"
#include "functions_AECM.hpp"
#include "functions_OSL.hpp"

#include "../CppNumericalSolvers/include/cppoptlib/meta.h"
#include "../CppNumericalSolvers/include/cppoptlib/boundedproblem.h"
#include "../CppNumericalSolvers/include/cppoptlib/solver/lbfgsbsolver.h"

#include <ctime>
#include <iomanip>









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
	r_scale = 22.0;
	r.array() /= r_scale;
	sigma.array() /= r_scale;
	Debug0("sigma: " << sigma.transpose());
	Debug2("Preprocessing done");
	
	
	

	Vector_eig TE_example((Vector_eig(3) << 0.03, 0.08, 0.01).finished());
	Vector_eig TR_example((Vector_eig(3) << 3, 3, 0.035).finished());
	
	double TE_scale = 2.01/TE_example.minCoeff();		// 1.01/0.03
	double TR_scale = 2.01/TR_example.minCoeff();		// 1.01/1.00
	Debug0("r_scale: " << r_scale);
	Debug0("TE scale: " << TE_scale);
	Debug0("TR scale: " << TR_scale);
	TE_example *= TE_scale;
	TR_example *= TR_scale;
	//TE_scale, TR_scale are needed for determining the bounds
	
	
	Vector_eig lb(3), ub(3);
	lb << 1e-6, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	for(int i = 1; i < 3; ++i){
		if(lb[i] < 1.0e-8){
			lb[i] = 1.0e-8;
		}
	}
	Debug0("lb:" << lb.transpose());
	Debug0("ub:" << ub.transpose());
	double W1_init = exp(-1/(2.0*TR_scale));		// exp(-1/(2.0*1.01))
	double W2_init = exp(-1/(0.1*TE_scale));		// exp(-1/(0.1*1.01/0.03))
	
	
	
	
	// Divide into train and test:
	std::vector<int> train_ind{0, 1, 2};
	std::vector<int> test_ind{0, 1, 2};
	
	
	// Least Sq:
	Matrix_eig_row W_orig = Init_val(r, TE_example, TR_example, our_dim, 
	                             r_scale, TE_scale, TR_scale, W1_init, W2_init, 1, will_write);
	Debug1("W_orig after noiseless LS: ");
	show_head(W_orig);
	std::cout << std::flush;
	
	
	// Write to a file: 
	std::ofstream file_LS;
	file_LS.open ("result/W_orig_29_brainweb.txt");
	for(int i = 0; i < W_orig.rows(); ++i){
		file_LS << W_orig.row(i) << "\n";
	}
	file_LS.close();
	
	
	
	
	
	
	
	
	
	
	
	
	// Generate new r matrix from W_orig: 
	
	// Change W_orig matrix due to TE_scale and TR_scale
	int n = W_orig.rows();
	for(int i = 0; i < n; ++i){
		// W_orig(i, 0) = W_orig(i, 0) * r_scale;
		W_orig(i, 1) = pow(W_orig(i, 1), 1./TR_scale);
		W_orig(i, 2) = pow(W_orig(i, 2), 1./TE_scale);
	}
	
	
	
	
	// New values: 
	// TE_example.clear(); TR_example.clear();	lb.clear(); ub.clear(); sigma.clear();
	
	TE_example.resize(12); TR_example.resize(12); sigma.resize(12);
	
	
	TE_example << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1;
	TR_example << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3;
	sigma << 1.99146, 1.81265, 1.82837, 2.30221, 1.63414, 1.71876, 3.13695, 1.77141, 1.55651, 2.72191, 1.63068, 1.4359;
	sigma.array() *= 10;
	Debug0("sigma: " << sigma.transpose());
	Matrix_eig generated_r = Gen_r(W_orig, TE_example, TR_example, sigma);
	
	
	TE_scale = 2.01/TE_example.minCoeff();		// 1.01/0.03
	TR_scale = 2.01/TR_example.minCoeff();		// 1.01/1.00
	Debug0("r_scale: " << r_scale);
	Debug0("TE scale: " << TE_scale);
	Debug0("TR scale: " << TR_scale);
	TE_example *= TE_scale;
	TR_example *= TR_scale;
	//TE_scale, TR_scale are needed for determining the bounds
	
	
	lb(0) = 1e-6; lb(1) = exp(-1/(0.01*TR_scale)); lb(2) = exp(-1/(0.001*TE_scale));
	ub(0) = 450.0; ub(1) = exp(-1/(4.0*TR_scale)); ub(2) = exp(-1/(0.2*TE_scale));
	for(int i = 1; i < 3; ++i){
		if(lb[i] < 1.0e-8){
			lb[i] = 1.0e-8;
		}
	}
	Debug0("lower bound:" << lb.transpose());
	Debug0("upper bound:" << ub.transpose());
	W1_init = exp(-1/(2.0*TR_scale));		// exp(-1/(2.0*1.01))
	W2_init = exp(-1/(0.1*TE_scale));		// exp(-1/(0.1*1.01/0.03))
	
		
	int m_total = 12;
	std::vector<int> whole_ind = {};
	for(int i = 0; i < m_total; ++i)
		whole_ind.push_back(i);
	
	train_ind = {0, 8, 9};
	test_ind.clear();
	test_ind = whole_ind;
	/*
	std::set_difference(whole_ind.begin(), whole_ind.end(), 
						train_ind.begin(), train_ind.end(),
                        std::inserter(test_ind, test_ind.begin()));
	*/
	
	std::cout << "\ntrain: \t";
	for(int n1 : train_ind) {
		std::cout << n1 << ' ';
	}
	std::cout << "\ntest: \t";
	for(int n1 : test_ind) {
		std::cout << n1 << ' ';
	}
	std::cout << "\n" << std::flush;
	
	
	
	Matrix_eig_row train(generated_r.rows(), train_ind.size());
	Vector_eig TE_train(train_ind.size()), TR_train(train_ind.size()), sigma_train(train_ind.size());
	short our_dim_train[8];
	for(int i = 0; i < train_ind.size(); ++i) {
		train.col(i) = generated_r.col(train_ind[i]);
		TE_train[i] = TE_example(train_ind[i]);
		TR_train[i] = TR_example(train_ind[i]);
		sigma_train[i] = sigma(train_ind[i]);
	}
	for(int i = 0; i < 8; ++i){
		our_dim_train[i] = our_dim[i];
	}
	our_dim_train[4] = (short)train_ind.size();		//our_dim[0] = 3 or 4
	
	Matrix_eig_row test(generated_r.rows(), test_ind.size());
	Vector_eig TE_test(test_ind.size()), TR_test(test_ind.size()), sigma_test(test_ind.size());
	for(int i = 0; i < test_ind.size(); ++i){
		test.col(i) = generated_r.col(test_ind[i]);
		TE_test[i] = TE_example(test_ind[i]);
		TR_test[i] = TR_example(test_ind[i]);
		sigma_test[i] = sigma(test_ind[i]);
	}
	
	Matrix_eig perf_1, perf_2, perf_3, perf_4;
	
	std::ofstream file_performance;
	file_performance.open ("result/Performances_29.txt");
	
	
	
	
	
	
	// W_orig
	for(int i = 0; i < n; ++i){
		// W_orig(i, 0) = W_orig(i, 0) * r_scale;
		W_orig(i, 1) = pow(W_orig(i, 1), TR_scale);
		W_orig(i, 2) = pow(W_orig(i, 2), TE_scale);
	}
	Debug1("W_orig transformed after noiseless LS: ");
	show_head(W_orig);
	std::cout << std::flush;
	
	
	
	
	// LS on the noise-contaminated image: 
	int do_least_sq = 1;
	Matrix_eig_row W_init = Init_val(train, TE_train, TR_train, our_dim_train, 
	                             r_scale, TE_scale, TR_scale, W1_init, W2_init, do_least_sq, will_write);
	Debug1("W_init after LS: ");
	show_head(W_init);
	std::cout << std::flush;
	
	
	Matrix_eig_row W_LS = W_init;
	
	
	perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
	Debug0("Avg perfs MPLE: " << perf_1.mean() << ", " << perf_2.mean() << ", "
						 << perf_3.mean() << ", " << perf_4.mean());
	std::cout << "Performances over images Penalized: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images Penalized: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images Penalized: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images Penalized: " << perf_4.transpose() << "\n\n\n" << std::endl;
	
	file_performance << "Performances over images Penalized: \t" << perf_1.transpose() << "\n";
	file_performance << "Performances over images Penalized: \t" << perf_2.transpose() << "\n";
	file_performance << "Performances over images Penalized: \t" << perf_3.transpose() << "\n";
	file_performance << "Performances over images Penalized: \t" << perf_4.transpose() << "\n\n\n";
	
	
	
	
	
	
	
	
	
	MRF_param MRF_obj_1(our_dim_train[1], our_dim_train[2], our_dim_train[3]);
	Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
	Vector_eig beta_init = 1.0*Vector_eig::Ones(3);
	
	
	
	// OSL: 
	OSL_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
	          our_dim_train[1], our_dim_train[2], our_dim_train[3], r_scale, TE_scale, TR_scale, MRF_obj_1, 
	          50, 1, 0.1, 1e-5, 1);
	//change
	Debug1("W - OSL");
	show_head(W_init);
	
	Matrix_eig_row W_OSL = W_init;
	
	// Write to a file: 
	std::ofstream file_OSL;
	file_OSL.open ("result/W_OSL_29_brainweb.txt");
	for(int i = 0; i < W_init.rows(); ++i){
		file_OSL << W_init.row(i) << "\n";
	}
	file_OSL.close();
	
		
	perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
	Debug0("Avg perfs MPLE: " << perf_1.mean() << ", " << perf_2.mean() << ", "
						 << perf_3.mean() << ", " << perf_4.mean());
	std::cout << "Performances over images Penalized: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images Penalized: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images Penalized: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images Penalized: " << perf_4.transpose() << "\n\n\n" << std::endl;
	
	file_performance << "Performances over images Penalized: \t" << perf_1.transpose() << "\n";
	file_performance << "Performances over images Penalized: \t" << perf_2.transpose() << "\n";
	file_performance << "Performances over images Penalized: \t" << perf_3.transpose() << "\n";
	file_performance << "Performances over images Penalized: \t" << perf_4.transpose() << "\n\n\n";
	
	
	
	
	
	
	
	
	// AECM: 
	
	AECM_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
	          our_dim_train[1], our_dim_train[2], our_dim_train[3], r_scale, TE_scale, TR_scale, MRF_obj_1, 
	          50, 1, 0.1, 1e-5, 1);
	//change
	Debug1("W - AECM");
	show_head(W_init);
	
	Matrix_eig_row W_AECM = W_init;
	
	// Write to a file: 
	std::ofstream file_final;
	file_final.open ("result/W_final_29_brainweb.txt");
	for(int i = 0; i < W_init.rows(); ++i){
		file_final << W_init.row(i) << "\n";
	}
	file_final.close();
	
		
	perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
	Debug0("Avg perfs MPLE: " << perf_1.mean() << ", " << perf_2.mean() << ", "
						 << perf_3.mean() << ", " << perf_4.mean());
	std::cout << "Performances over images Penalized: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images Penalized: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images Penalized: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images Penalized: " << perf_4.transpose() << "\n\n\n" << std::endl;
	
	file_performance << "Performances over images Penalized: \t" << perf_1.transpose() << "\n";
	file_performance << "Performances over images Penalized: \t" << perf_2.transpose() << "\n";
	file_performance << "Performances over images Penalized: \t" << perf_3.transpose() << "\n";
	file_performance << "Performances over images Penalized: \t" << perf_4.transpose() << "\n\n\n";
	file_performance.close();
	
	
	
	
	
	
	std::time_t t2 = std::time(nullptr);
	std::tm tm2 = *std::localtime(&t2);
	std::cout << "Current time: " << std::put_time(&tm2, "%c %Z") << '\n';

	return 0;
}



