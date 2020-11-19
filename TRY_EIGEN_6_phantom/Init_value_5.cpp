/**
*




https://slideplayer.com/slide/5272436/
alpha = 90 degree would give this equation (http://mriquestions.com/spoiled-gre-parameters.html)
(But also see parameter selection table)
(Also, I have no idea are these related or not: http://www.mriquestions.com/bloch-equations.html)
http://www.mriquestions.com/image-contrast-trte.html
http://www.mriquestions.com/complete-list-of-questions.html





g++ Init_value_5.cpp -o test_LS -I /usr/include/eigen3 -O3 -lgsl -lgslcblas -lm

./test_LS ../Read_Data/new_phantom.nii Dummy_sd_3D.txt 0

./test_LS ../Read_Data/small_phantom.nii Dummy_sd.txt 0


*/






#include "scheme_new_numerical.hpp"
#include "Read_files_2.hpp"



#include "../CppNumericalSolvers/include/cppoptlib/meta.h"
#include "../CppNumericalSolvers/include/cppoptlib/boundedproblem.h"
#include "../CppNumericalSolvers/include/cppoptlib/solver/lbfgsbsolver.h"










template<typename T>
class Least_Sq_est : public cppoptlib::BoundedProblem<T> {
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	using TMatrix = typename cppoptlib::BoundedProblem<T>::THessian;
	typedef Matrix_eig_row TMatrix_row;
	
	TMatrix_row r;
	TVector r2;


  public:
	Least_Sq_est(const TVector y_) : 
		cppoptlib::BoundedProblem<T>(y_.size()), r2(y_){}
	


	TVector TE, TR, lb, ub;
	int i;
	
	// Track the best:
	Eigen::VectorXd current_best_param;
	double current_best_val = 1.0e+15;


	// Objective function, to be minimized:
	T value(const TVector &x) {
		TVector v_new = Bloch_vec(x, TE, TR);
		// if(i == 0)
		if(i == 2){				// corresponding to 3 as in R
			DebugLS("x: " << x.transpose() << "; value: " << (r.row(i).transpose() - v_new).squaredNorm());
		}
		double fx = (r.row(i).transpose() - v_new).squaredNorm();
		
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
		TVector v_new = Bloch_vec(x, TE, TR);
		grad << 0,0,0;
		int m = TR.size();
		
		
		for(int j = 0; j < m; ++j){
		
			grad[0] -= 2*(r(i, j) - v_new(j)) * 
						std::exp(TE(j)*std::log(x(2))) * 
						(1-std::exp(TR(j)*std::log(x(1))));
			
			grad[1] -= 2*(v_new(j) - r(i, j)) * 
						x(0)*TR(j) * std::exp(TE(j)*log(x(2))) * 
						std::exp((TR(j)-1)*std::log(x(1)));
			
			grad[2] -= 2*(r(i, j) - v_new(j)) * 
						x(0)*TE(j) * std::exp((TE(j)-1)*log(x(2))) * 
						(1-std::exp(TR(j)*std::log(x(1))));
		}
		//if(i == 0)
		if(i == 2)		// corresponding to 3 as in R
			DebugLS("grad: " << grad.transpose() << "\n" );
	}
};





/*
* Least square solution
* 	Input:  
*			W, TE_example, TR_example, r, r_scale, TE_scale, TR_scale
* 			
* 			W is changed
*/
void least_sq_solve(Matrix_eig_row &W, 
					const Eigen::VectorXd &TE_example, const Eigen::VectorXd &TR_example, 
                    const Matrix_eig_row &r, double r_scale, double TE_scale, double TR_scale){


	Debug0("Doing Least Square Estimate!");
	auto time_1_lsq = std::chrono::high_resolution_clock::now();
	
	Eigen::VectorXd r2 = Eigen::VectorXd::Ones(3) * 0.5;
	Least_Sq_est<double> f(r2);
	Eigen::VectorXd x(3), lb(3), ub(3); 
	
	//Bounds of rho, W1, W2:
	lb << 0.0001, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	for(int i = 1; i < 3; ++i){
		if(lb[i]<1.0e-8){
			lb[i] = 1.0e-8;
		}
	}
	Debug1("lb inside LS: " << lb.transpose());
	Debug1("ub inside LS: " << ub.transpose());
	
	
	f.r.noalias() = r;	f.TE.noalias() = TE_example;	f.TR.noalias() = TR_example;
	f.setLowerBound(lb);	f.setUpperBound(ub);		f.lb.noalias() = lb; 	f.ub.noalias() = ub;
	
	double old_val = 0.0, fx;
	int n = r.rows(), bad_count_o = 0, bad_count_o_2 = 0, bad_bound_1 = 0, bad_bound_2 = 0, nan_count = 0;



	
	
	// Declaring the solver
	cppoptlib::LbfgsbSolver<Least_Sq_est<double>> solver;
	cppoptlib::Criteria<double> crit_LS = cppoptlib::Criteria<double>::defaults();
	crit_LS.iterations = 100;
	solver.setStopCriteria(crit_LS);
	
	
	
	
	// Loop of 
	for(int i = 0; i < n; ++i){
	
		if(i==100000 || i==200000 || i==300000 || i==400000 || i==500000 || i==600000 || i==700000 || i==800000 || i==900000 ){
			Debug1("i: "<< i);
		}
		
		// Track the best:
		// double current_best_val = 1.0e+15;
		f.current_best_val = 1.0e+15;		// This was not updated previously   Nov 8, 2.36 am
		
		
		f.i = i;
		x = W.row(i);
		
		//Check Bounds:
		/*
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
		*/
		
		
	
		//Print initial values:
		DebugLS ("value of i: " << i << "\t x at first: " << x.transpose());
		old_val = f.value(x);
		DebugLS ("f(x) at first: \n" << old_val) ;
		
	
		//Solve:
		solver.minimize(f, x);
		fx = f(x);
		
		// Track the best:
		x = f.current_best_param;
		fx = f.current_best_val;
		DebugLS("f(param_new) in argmin: " << fx << "\t x:" << x.transpose());
		
		//if(i == 0)
		if(i == 2){			// corresponding to 3 as in R
			DebugLS("argmin: " << x.transpose() << ";\tf(x) in argmin: " << f(x)) ;
			DebugLS("Solver status: " << solver.status() );	//Guess: bad reports: under constraints => grad is not ~0 
			DebugLS("Final criteria values: " << "\n" << solver.criteria());
			// DebugLS("f(param_new) in argmin: " << fx << "\t old val:" << old_val);
		}
		
		
		
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
	auto duration_lsq = std::chrono::duration_cast<std::chrono::seconds>(time_2_lsq - time_1_lsq);
	Debug1("Time taken total Least Square: " << duration_lsq.count() << " seconds\n");
}






/*
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





/*
* Creates the initial matrix W
* if do_least_sq is 1, 
* 	it gives the least square solution.
*/
Matrix_eig_row Init_val(const Matrix_eig_row &r, 
                        const Eigen::VectorXd &TE_example, const Eigen::VectorXd &TR_example, 
                        short our_dim[8], 
                        double r_scale, double TE_scale, double TR_scale, 
                        double W_1_init = exp(-1/2.0), double W_2_init = exp(-1/0.1), 
                        int do_least_sq = 1, char will_write = 0){


	//Primary Initial value for test//
	int n = our_dim[1]*our_dim[2]*our_dim[3];
	Matrix_eig_row W = Matrix_eig_row::Ones(n, 3);
	//show_dim(W);
	//show_dim(r);
	
	
	// Ad hoc initial values:
	W.col(0) = r.rowwise().mean().transpose();
	W.col(1) *= W_1_init;
	W.col(2) *= W_2_init;
	for(int i = 0; i < r.rows(); ++i){
		//if(W(i, 0)>450.0/r_scale){
		//	W(i, 0) = 425.0/r_scale;
		//}
		if(W(i, 0)>450.0){
			W(i, 0) = 425.0;
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






int main(int argc, char * argv[]) {

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
	
	
	
	
	
	
	


	//Vector_eig TE_example((Vector_eig(12) << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1).finished());
	//Vector_eig TR_example((Vector_eig(12) << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3).finished());




	Vector_eig TE_example((Vector_eig(18) << 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10).finished());
	Vector_eig TR_example((Vector_eig(18) << 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3).finished());
	// 1.01 -> 2.01
	double TE_scale = 2.01/TE_example.minCoeff();		// 1.01/0.03
	double TR_scale = 2.01/TR_example.minCoeff();		// 1.01/1.00
	// Debug0("r_scale: " << r_scale);
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
	
	//std::vector<int> train_ind{0, 1, 2};
	//std::vector<int> test_ind{3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
	
	std::vector<int> train_ind{0, 6, 13};
	std::vector<int> test_ind{1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17};
	// Also, this creates 450 in the first edge - but not with non-penalized case - check
	// Somehow only 0, 1, 2 in trainset creates Nan's. We have to look.
	
	
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
	
	std::ofstream file_performance;
	file_performance.open ("result/Performances_LS.txt");


	
	
	// Temp results: Performance on the Init W: 
	
	
	// Least Sq:
	// Change 
	int do_least_sq = 1;	// 0 Subrata -- least sq have better initial likelihood-but stucks and gives nan in some value
	Matrix_eig_row W_init = Init_val(train, TE_train, TR_train, our_dim_train, 
	                             r_scale, TE_scale, TR_scale, W1_init, W2_init, do_least_sq, will_write);
	Debug1("W initial done");
	check_nan(W_init, "W matrix init, nan: \n");
	// int nan_count_1st = check_nan_W(W_init, W_1st);		// Change as there is no W_1st
	// Debug0("NAN count at first:" << nan_count_1st);
	Debug1("W_init after LS: ");
	show_head(W_init);
	std::cout << std::flush;
	
	
	
}



















/*
int main(int argc, char * argv[]){
	
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
	// r.array() /= r_scale;
	// sigma.array() /= r_scale;
	
	Debug0("sigma: " << sigma.transpose());
	Debug2("Preprocessing done");
	
	
	
	
	
	
	


	//Vector_eig TE_example((Vector_eig(12) << 0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1).finished());
	//Vector_eig TR_example((Vector_eig(12) << 0.6, 0.6, 0.6, 1, 1, 1, 2, 2, 2, 3, 3, 3).finished());




	Vector_eig TE_example((Vector_eig(18) << 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10).finished());
	Vector_eig TR_example((Vector_eig(18) << 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3).finished());
	// 1.01 -> 2.01
	double TE_scale = 2.01/TE_example.minCoeff();		// 1.01/0.03
	double TR_scale = 2.01/TR_example.minCoeff();		// 1.01/1.00
	// Debug0("r_scale: " << r_scale);
	Debug0("TE scale: " << TE_scale);
	Debug0("TR scale: " << TR_scale);
	TE_example *= TE_scale;
	TR_example *= TR_scale;
	//TE_scale, TR_scale are needed for determining the bounds
	
	
	
	Vector_eig lb(3), ub(3);
	//lb << 0.0001/r_scale, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	//ub << 450.0/r_scale, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	
	lb << 0.0001, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale));
	ub << 450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale));
	
	Debug0("lb:" << lb.transpose());
	Debug0("ub:" << ub.transpose());
	
	
	
	double W1_init = exp(-1/(2.0*TR_scale));		// exp(-1/(2.0*1.01))
	double W2_init = exp(-1/(0.1*TE_scale));		// exp(-1/(0.1*1.01/0.03))
	




	
	// Divide into train and test:
	
	//std::vector<int> train_ind{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	//std::vector<int> test_ind{10, 11};

	//std::vector<int> train_ind{0, 9, 11};
	//std::vector<int> test_ind{1, 2, 3, 4, 5, 7, 8, 10};
        
        
	//std::vector<int> train_ind{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	//std::vector<int> test_ind{16, 17};
	
	std::vector<int> train_ind{0, 1, 2};
	std::vector<int> test_ind{3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
	// Also, this creates 450 in the first edge - but not with non-penalized case - check
	// Somehow only 0, 1, 2 in trainset creates Nan's. We have to look.
	
	
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
	
	std::ofstream file_performance;
	file_performance.open ("result/Performances_17_new.txt");


	
	
	// Temp results: Performance on the Init W: 
	
	
	// Least Sq:
	// Change 
	int do_least_sq = 1;	// 0 Subrata -- least sq have better initial likelihood-but stucks and gives nan in some value
	Matrix_eig_row W_init = Init_val(train, TE_train, TR_train, our_dim_train, 
	                             r_scale, TE_scale, TR_scale, W1_init, W2_init, do_least_sq, will_write);
	Debug1("W initial done");
	check_nan(W_init, "W matrix init, nan: \n");
	// int nan_count_1st = check_nan_W(W_init, W_1st);		// Change as there is no W_1st
	// Debug0("NAN count at first:" << nan_count_1st);
	Debug1("W_init after LS: ");
	show_head(W_init);
	std::cout << std::flush;
	
	
	
	
	//https://eigen.tuxfamily.org/dox-devel/group__TutorialSlicingIndexing.html
	// However, 'all' is still not in the current stable version.
	// https://stackoverflow.com/questions/58699992/why-eigen-does-not-resolve-built-in-symbols-all-last-seq-etc



	
	
	
	//Primary Initial value for test//
	int n = our_dim_train[1]*our_dim_train[2]*our_dim_train[3];
	Matrix_eig_row W = Matrix_eig_row::Ones(n, 3);
	
	// Ad hoc initial values:
	// W.col(0) = r.rowwise().mean().transpose();
	W.col(0) = train.rowwise().mean().transpose();	// So, nan or not - depends on inital value???
	W.col(1) *= W1_init;
	W.col(2) *= W2_init;
	for(int i = 0; i < train.rows(); ++i){
		if(W(i, 0)>450){
			W(i, 0) = 425;
		}
	}


	show_head(W);
	Debug0("W.row(55643): " << W.row(55643));	
	Debug0("train.row(55643): " << train.row(55643));
	least_sq_solve(W, TE_train, TR_train, train, r_scale, TE_scale, TR_scale);
		
	
	
	
	
	
	
	// Least Sq:
	int do_least_sq = 1;	// 0 Subrata -- least sq have better initial likelihood-but stucks and gives nan in some value
	Matrix_eig_row W_init = Init_val(train, TE_train, TR_train, our_dim_train, 
	                             r_scale, TE_scale, TR_scale, W1_init, W2_init, do_least_sq, will_write);
	Debug1("W initial done");
	check_nan(W_init, "W matrix nan initially: \n");
	// int nan_count_1st = check_nan_W(W_init, W_1st);		// Change as there is no W_1st
	// Debug0("NAN count at first:" << nan_count_1st);
	Debug1("W_init after LS: ");
	show_head(W_init);
	std::cout << std::flush;
	
	return 0;
}

*/
