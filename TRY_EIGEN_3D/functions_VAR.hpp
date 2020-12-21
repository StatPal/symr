


/* _VAR_HEADER_ */
#ifndef _VAR_HEADER_
#define _VAR_HEADER_




#include "scheme_new_numerical.hpp"
#include "Init_value_6_numerical.hpp"
#include "functions_AECM.hpp"

#include <ctime>
#include <iomanip>










/*
* Hessian matrix iterative solution:
* v_grad ' Hessian_mat_without_MRF v_grad is to be calculated 
* j-th image's variance(n = n_x*n_y*n_z) is to be calculated.
*/
Matrix_eig Var_est_test_mat(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                            const Vector_eig &TE_train, const Vector_eig &TR_train,
                            const Vector_eig &sigma_train, const Matrix_eig_row &train, 
                            int n_x, int n_y, int n_z, MRF_param &MRF_obj,
                            const Vector_eig &TE_test, const Vector_eig &TR_test, 
                            const Vector_eig &sigma_test, const Matrix_eig_row &test){

	auto hess_1 = std::chrono::high_resolution_clock::now();
	int n = W.rows();
	Matrix_eig Var_est(n, TE_test.size());
	
	SpMat A = Hessian_mat(W, Psi_inv, beta, TE_train, TR_train, sigma_train, train, n_x, n_y, n_z, MRF_obj, 1);	
	//SpMat A = Hessian_mat(W, Psi_inv, beta, TE_train, TR_train, sigma_train, train, n_x, n_y, n_z, MRF_obj, 0);
	assert(A.rows() == 3*n);
	save_sparse(A, "Hessian_Matrix_26.csv", 1);
	
	// Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
	// Debug0("Diagonal");
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
	Debug0("IncompleteCholesky");
	cg.compute(A);
	
	
	
	
	// Vector_eig tmp_soln(n);
	// SpVec b(n);
	// Wait, this is a BUG - how did this went unnoticed!!!
	Vector_eig tmp_soln(3*n);
	SpVec b(3*n);
	// Vector_eig b = Vector_eig::Zero(3*n);
	
	
	// First i or j, which would be better? - check.
	std::cout << std::endl << std::flush;
	Debug1("Hessian loop starts!");
	std::cout << std::endl << std::flush;
	for(int j = 0; j < TE_test.size(); ++j){
		for(int i = 0; i < n; ++i){
		
			//if(i==100000 || i==300000 || i==500000 || i==700000 || i==900000 ){
			if( i % 100 == 0){
				std::cout << std::endl;
				Debug1("Info i: "<< i << ", j: " << j);
			}
			
			// b = v_grad(W, Psi_inv, beta, TE_test, TR_test, sigma_test, test, n_x, n_y, n_z, i, j);
			v_grad(W, Psi_inv, beta, TE_test, TR_test, sigma_test, test, n_x, n_y, n_z, i, j, b);
			
			assert(b.size() == A.cols());
			
			// Thought about the j
			// Should it be ind [j] or something like that?
			// I don't think so!
			tmp_soln = cg.solve(b);
			if( i % 100 == 0)
				std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
			Var_est(i, j) = b.dot(tmp_soln);
			
			// b(3*i) = 0.0; b(3*i+1) = 0.0; b(3*i+2) = 0.0;
			
		}
	}
	auto hess_2 = std::chrono::high_resolution_clock::now();
	auto hess_duration = std::chrono::duration_cast<std::chrono::seconds>(hess_2 - hess_1);
	Debug1("Time taken for Info matrix using Hessian: " << hess_duration.count() << " seconds\n");
	
	return Var_est;
}







/*
* Parametric Bootstrap
* 'const' are removed as OSL_optim inside needs non-const cases
* -- No, just create another set of small matrices. Otherwise in main, Psi and beta would be changed
* wait, there are no terain - test case?
*/


/*
* Parametric Bootstrap
*/
Matrix_eig para_boot_test_mat(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                              const Vector_eig &TE_train, const Vector_eig &TR_train, 
                              const Vector_eig &sigma_train, const Matrix_eig_row &r,
                              int n_x, int n_y, int n_z, 
                              double r_scale, double TE_scale, double TR_scale, MRF_param &MRF_obj, 
                              const Vector_eig &TE_test, const Vector_eig &TR_test, 
                              const Vector_eig &sigma_test, const Matrix_eig_row &test,
                              int B = 15, int EM_iter = 10, double abs_diff = 1.0e-4, double rel_diff = 1e-3){

	auto boot_1 = std::chrono::high_resolution_clock::now();
	Debug1("Parametric Bootstrap starts!!");
	int n = W.rows();
	int m_test = TE_test.size();
	int m_train = TE_train.size();
	Matrix_eig Var_est(n, m_test);
	Matrix_eig_row W_init = W;
	Matrix3d_eig Psi_inv_init = Psi_inv;
	Vector_eig beta_init = beta;
	
	// All estimated parametrs:
	Matrix_eig generated_r(n, m_train);		// train columns only
	
	// added cases: 
	Matrix_eig tmp_mat = Matrix_eig::Zero(n, m_test);
	Matrix_eig sum_mat = Matrix_eig::Zero(n, m_test);
	Matrix_eig sum_sq_mat = Matrix_eig::Zero(n, m_test);
	
	for(int b = 0; b < B; ++b){
		if( b % 50 == 0){
			Debug0("bootstrap sample " << b);
		}
		// check_nan(W, "W matrix in boot, nan: \n");

		
		//Generate an image matrix:
		// generated_r = Gen_r(W, TE_test, TR_test, sigma_test); // Sorry, this is a mistake
		generated_r = Gen_r(W, TE_train, TR_train, sigma_train);
		// W_init.noalias() = W;		// Not needed? - numerical stabilty?
		OSL_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, generated_r, 
							n_x, n_y, n_z, r_scale, TE_scale, TR_scale, MRF_obj, EM_iter, 1, abs_diff, rel_diff, 0);
		tmp_mat = v_mat(W_init, TE_test, TR_test);
		sum_mat += tmp_mat;
		sum_sq_mat += tmp_mat.array().square().matrix();
	}
	
	sum_mat /= B;
	sum_sq_mat /= B;
	sum_sq_mat -= sum_mat.array().square().matrix();
	
	auto boot_2 = std::chrono::high_resolution_clock::now();
	auto boot_duration = std::chrono::duration_cast<std::chrono::seconds>(boot_2 - boot_1);
	Debug1("Time taken for Info matrix using Parametric Bootstrap: " << boot_duration.count() << " seconds\n");
	
	return (sum_sq_mat);
}





#endif /* !_VAR_HEADER_ */
