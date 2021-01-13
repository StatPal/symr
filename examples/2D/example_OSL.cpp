/**
* 
* OSL EM algorithm
* Psi and beta are updated at every loop

* E step done seperately : 




* To compile:

g++ example_OSL.cpp -o example_OSL -I /usr/include/eigen3 -O3 -lgsl -lgslcblas -lm -fopenmp -DEIGEN_DONT_PARALLELIZE


./example_OSL ../Read_Data/new_phantom.nii Dummy_sd.txt 0

./example_OSL ../../data/new_phantom.nii Dummy_sd.txt 0

./example_OSL ../Read_Data/small_phantom.nii Dummy_sd.txt 0

nohup ./example_OSL ../Read_Data/new_phantom.nii Dummy_sd.txt 0 > example_OSL.out & 




* 
*/



#include "../../include/2D/read_files.hpp"
#include "../../include/2D/functions_gen.hpp"
#include "../../include/2D/functions_LS_and_init_value.hpp"
#include "../../include/2D/functions_OSL.hpp"

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
	r_scale = 10.0;
	r.array() /= r_scale;
	sigma.array() /= r_scale;
	
	Debug0("sigma: " << sigma.transpose());
	Debug2("Preprocessing done");
	
	
	
	
	
	
	
	
	
	Vector_eig TE_example((Vector_eig(18) << 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10).finished());
	Vector_eig TR_example((Vector_eig(18) << 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3).finished());
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
	




	
	// Divide into train and test:
//	std::vector<int> train_ind{0, 6, 13};
//	std::vector<int> test_ind{1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17};
	
	int m_total = 18;
	std::vector<int> whole_ind = {};
	for(int i = 0; i < m_total; ++i)
		whole_ind.push_back(i);
	
	std::vector<int> train_ind{0, 6, 13};
	std::vector<int> test_ind{};
	
	test_ind = whole_ind;
	/*
	std::set_difference(whole_ind.begin(), whole_ind.end(), 
						train_ind.begin(), train_ind.end(),
                        std::inserter(test_ind, test_ind.begin()));
	*/
	
	
	
	Matrix_eig_row train(r.rows(), train_ind.size());
	Vector_eig TE_train(train_ind.size()), TR_train(train_ind.size()), sigma_train(train_ind.size());
	short our_dim_train[8];
	for(int i = 0; i < (int)train_ind.size(); ++i) {
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
	for(int i = 0; i < (int)test_ind.size(); ++i){
		test.col(i) = r.col(test_ind[i]);
		TE_test[i] = TE_example(test_ind[i]);
		TR_test[i] = TR_example(test_ind[i]);
		sigma_test[i] = sigma(test_ind[i]);
	}
	
	Matrix_eig perf_1, perf_2, perf_3, perf_4;
	
	std::ofstream file_performance;
	file_performance.open ("result/Performances_26.txt");


	
	
	
	
	
	
	// Least Sq:
	// Change 
	int do_least_sq = 1;
	Matrix_eig_row W_init = Init_val(train, TE_train, TR_train, our_dim_train, 
	                             r_scale, TE_scale, TR_scale, W1_init, W2_init, do_least_sq, will_write);
	Debug1("W_init after LS: ");
	show_head(W_init);
	std::cout << std::flush;
	
	
	// Write to a file: 
	std::ofstream file_LS;
	file_LS.open ("result/W_LS_26.txt");
	for(int i = 0; i < W_init.rows(); ++i){
		file_LS << W_init.row(i) << "\n";
	}
	file_LS.close();
	
	
	// Save the estimated values: 
	Vector_eig v_new = Vector_eig::Zero(TE_test.size());
	std::ofstream file_predicted;
	file_predicted.open ("result/v_predicted_LS_26.txt");
	for(int i = 0; i < W_init.rows(); ++i){
		Bloch_vec(W_init.row(i), TE_test, TR_test, v_new);
		file_predicted << v_new.transpose() << "\n";
	}
	file_predicted.close();
	
	
	
	perf_1 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2);
	Debug0("Avg perfs LS: " << perf_1.mean() << ", " << perf_2.mean() << ", "
						 << perf_3.mean() << ", " << perf_4.mean());
	std::cout << "Performances over images LS: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images LS: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images LS: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images LS: " << perf_4.transpose() << "\n\n\n" << std::endl;
	
	file_performance << "Performances over images LS: \t" << perf_1.transpose() << "\n";
	file_performance << "Performances over images LS: \t" << perf_2.transpose() << "\n";
	file_performance << "Performances over images LS: \t" << perf_3.transpose() << "\n";
	file_performance << "Performances over images LS: \t" << perf_4.transpose() << "\n\n\n";
	
	
	
	
	
	
	
	
	// Test:
	Matrix_eig_row W_LS = W_init;
	

	
	
	
	// Likelihood Based optimization:
	Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
	Vector_eig beta_init = 1.0*Vector_eig::Ones(3);
	
	MRF_param MRF_obj_1(our_dim_train[1], our_dim_train[2], our_dim_train[3]);
	
	
	
	
	// Non -penalized:
	
	OSL_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
	          our_dim_train[1], our_dim_train[2], our_dim_train[3], r_scale, TE_scale, TR_scale, MRF_obj_1, 
	          50, 0, 0.1, 1e-5, 1);
	//change
	
	// Write to a file: 
	std::ofstream file_Likeli;
	file_Likeli.open ("result/W_Likeli_26.txt");
	for(int i = 0; i < W_init.rows(); ++i){
		file_LS << W_init.row(i) << "\n";
	}
	file_Likeli.close();
	
	
	// Save the estimated values: 
	v_new = Vector_eig::Zero(TE_test.size());
	file_predicted.open ("result/v_predicted_Likeli_26.txt");
	for(int i = 0; i < W_init.rows(); ++i){
		Bloch_vec(W_init.row(i), TE_test, TR_test, v_new);
		file_predicted << v_new.transpose() << "\n";
	}
	file_predicted.close();
	
	
	Matrix_eig_row W_likeli = W_init;
	W_init.noalias() = W_LS;
	show_head(W_likeli);
	
	
	
	
	
	perf_1 = Performance_test(W_likeli, test, TE_test, TR_test, sigma_test, 1, 1);
	perf_2 = Performance_test(W_likeli, test, TE_test, TR_test, sigma_test, 3, 1);
	perf_3 = Performance_test(W_likeli, test, TE_test, TR_test, sigma_test, 1, 2);
	perf_4 = Performance_test(W_likeli, test, TE_test, TR_test, sigma_test, 3, 2);
	Debug0("Avg perfs MLE: " << perf_1.mean() << ", " << perf_2.mean() << ", "
						 << perf_3.mean() << ", " << perf_4.mean());
	std::cout << "Performances over images Likelihood: " << perf_1.transpose() << "\n";
	std::cout << "Performances over images Likelihood: " << perf_2.transpose() << "\n";
	std::cout << "Performances over images Likelihood: " << perf_3.transpose() << "\n";
	std::cout << "Performances over images Likelihood: " << perf_4.transpose() << "\n\n\n" << std::endl;
	
	
	file_performance << "Performances over images Likelihood: \t" << perf_1.transpose() << "\n";
	file_performance << "Performances over images Likelihood: \t" << perf_2.transpose() << "\n";
	file_performance << "Performances over images Likelihood: \t" << perf_3.transpose() << "\n";
	file_performance << "Performances over images Likelihood: \t" << perf_4.transpose() << "\n\n\n";
	
	
	
	
	
	
	
	// Penalised:
	
	OSL_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
	          our_dim_train[1], our_dim_train[2], our_dim_train[3], r_scale, TE_scale, TR_scale, MRF_obj_1, 
	          50, 1, 0.1, 1e-5, 1);
	//change
	Debug1("W - Penalized Likelihood");
	show_head(W_init);
	
	// Write to a file: 
	std::ofstream file_final;
	file_final.open ("result/W_final_26.txt");
	for(int i = 0; i < W_init.rows(); ++i){
		file_final << W_init.row(i) << "\n";
	}
	file_final.close();
	
	// Save the estimated values: 
	v_new = Vector_eig::Zero(TE_test.size());
	file_predicted.open ("result/v_predicted_OSL_26.txt");
	for(int i = 0; i < W_init.rows(); ++i){
		Bloch_vec(W_init.row(i), TE_test, TR_test, v_new);
		file_predicted << v_new.transpose() << "\n";
	}
	file_predicted.close();
	
	
	
	
	
	
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






