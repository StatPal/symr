/**
* 
Example variance for contrast vector :




* To compile :

g++ example_VAR_part.cpp -o example_VAR_part -I /usr/include/eigen3 -O3 -lgsl -lgslcblas -lm -fopenmp

g++ example_VAR_part.cpp -o example_VAR_part -I ~/program/eigen3 -O3 -lgsl -lgslcblas -lm -fopenmp -DEIGEN_DONT_PARALLELIZE


./example_VAR_part ../Read_Data/ZHRTS2.nii Dummy_sd_3D.txt 0

./example_VAR_part ../../data/ZHRTS2.nii Dummy_sd_3D.txt 0

./example_VAR_part ../Read_Data/small.nii Dummy_sd_3D.txt 0

nohup ./example_VAR_part ../Read_Data/ZHRTS2.nii Dummy_sd_3D.txt 0 > example_VAR.out & 

* 
*/




#include "../../include/3D/functions_gen.hpp"
#include "../../include/3D/read_files.hpp"
#include "../../include/3D/functions_LS_and_init_value.hpp"

#include "../../include/3D/functions_AECM.hpp"
#include "../../include/3D/functions_VAR.hpp"

#include <ctime>
#include <iomanip>










int main(int argc, char * argv[]) {

	
	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);
	std::cout << "Current time: " << std::put_time(&tm, "%c %Z") << "\n";
		
	if (argc != 5) {
		fprintf(stderr, "\nUsage: %s <data_file_name> <class_file_name> <SD_file_name> <will_write_to_a_file?> \n", argv[0]);
		exit(EXIT_FAILURE);
	}
	char *data_file, *sd_file, *class_file;
	data_file = argv[1]; 	class_file = argv[2];	sd_file = argv[3];
	char will_write = *(argv[4])-48;		// Converted from ascii
	short our_dim[8];
	
	
	
	
	// Reading the data: 
	Matrix_eig_row r = Preprocess_data(data_file, our_dim, will_write);
	Vector_eig sigma = read_sd(sd_file, our_dim[4]);
	
	Matrix_eig_row contrast = Read_nift1(class_file, our_dim, will_write);
	//Debug1("Contrast head");
	//show_head(contrast);
	
	
	// Scaled: r, sigma, ub would change.
	double r_scale = r.maxCoeff();
	r_scale = 1.0;
	r.array() /= r_scale;
	sigma.array() /= r_scale;
	
	Debug0("sigma: " << sigma.transpose());
	Debug2("Preprocessing done");
	
	
	
	Vector_eig TE_example((Vector_eig(12) << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1).finished());
	Vector_eig TR_example((Vector_eig(12) << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3).finished());
	
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
		if(lb[i] < 1.0e-8){
			lb[i] = 1.0e-8;
		}
	}
	Debug0("lower bound:" << lb.transpose());
	Debug0("upper bound:" << ub.transpose());
	double W1_init = exp(-1/(2.0*TR_scale));		// exp(-1/(2.0*1.01))
	double W2_init = exp(-1/(0.1*TE_scale));		// exp(-1/(0.1*1.01/0.03))
	




	
	
	
	
	// Divide into train and test:
	
	int m_total = 12;
	std::vector<int> whole_ind = {};
	for(int i = 0; i < m_total; ++i)
		whole_ind.push_back(i);
	
	std::vector<int> train_ind{0, 8, 9};
	std::vector<int> test_ind{};
	
	//test_ind = whole_ind;
	
	std::set_difference(whole_ind.begin(), whole_ind.end(), 
						train_ind.begin(), train_ind.end(),
                        std::inserter(test_ind, test_ind.begin()));
	
	
	std::cout << "\ntrain: \t";
	for(int n1 : train_ind) {
		std::cout << n1 << ' ';
	}
	std::cout << "\ntest: \t";
	for(int n1 : test_ind) {
		std::cout << n1 << ' ';
	}
	std::cout << "\n" << std::flush;
	
	
	
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
	
	Matrix_eig perf_1, perf_2, perf_3, perf_4;
	
	
	
	
	
	
	
	
	
	
	
	// Least Sq:
	// Change 
	int do_least_sq = 1;
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
	
	
	
	
	// Likelihood Based optimization:
	Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
	Vector_eig beta_init = 1.0*Vector_eig::Ones(3);
	MRF_param MRF_obj_1(our_dim_train[1], our_dim_train[2], our_dim_train[3]);
	// Penalised:
	AECM_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
	          our_dim_train[1], our_dim_train[2], our_dim_train[3], r_scale, TE_scale, TR_scale, MRF_obj_1, 
	          500, 1, 0.1, 1e-8, 1);
	//change
	Debug1("W - Penalized Likelihood");
	show_head(W_init);
	
	perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
	Debug0("Avg perfs MPLE: " << perf_1.mean() << ", " << perf_2.mean() << ", "
						 << perf_3.mean() << ", " << perf_4.mean());
	
	
	
	
	
	
	// The contrast: 
	int n = W_init.rows();
	
	SpVec contrast_2(n);
	for(int i = 0; i < n; ++i){
		if(contrast(i) == 3){
			contrast_2.insert(i) = 1.0;
		}
	}
	
	
	
	/* Variance estimation: */ 
	
	// Using Info matrix + delta method:	
	Vector_eig info_var_contrast = Var_est_test_mat_contrast(W_init, Psi_inv_init, beta_init, 
											 TE_train, TR_train, sigma_train,  
                                             train, our_dim_train[1], our_dim_train[2], our_dim_train[3], MRF_obj_1,
                                             TE_test, TR_test, sigma_test, test, contrast_2);
	
	
	
	// Using Bootstrap
	Vector_eig boot_var_contrast = para_boot_test_mat_contrast(W_init, Psi_inv_init, beta_init, 
											   TE_train, TR_train, sigma_train,  
                                               train, our_dim_train[1], our_dim_train[2], our_dim_train[3],
                                               r_scale, TE_scale, TR_scale, MRF_obj_1,
                                               TE_test, TR_test, sigma_test, test, contrast_2,
                                               200, 50, 1e-1, 1e-3);
	
	
	
	std::cout << "\n\nVariance from Information matrix:\n";
	Debug0("Info var contrast : " << info_var_contrast.transpose());
	
	
	std::cout << "\nVariance from parametric bootstrap:\n";
	Debug0("boot_var_contrast: " << boot_var_contrast.transpose());
	
	
	
	
	std::time_t t2 = std::time(nullptr);
	std::tm tm2 = *std::localtime(&t2);
	std::cout << "Current time: " << std::put_time(&tm2, "%c %Z") << '\n';

	return 0;
}


