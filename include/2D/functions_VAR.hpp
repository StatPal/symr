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


/* _VAR_HEADER_ */
#ifndef _VAR_HEADER_
#define _VAR_HEADER_




#include "../functions_gen.hpp"
#include "../functions_LS_and_init_value.hpp"
#include "../functions_AECM.hpp"

#include <ctime>
#include <iomanip>









/***************************************************
**************** Information Matrix ****************
****************************************************/


/**
* @brief	Derivative of \f$I_1/I_0\f$
* @param	x
* @return	\f$d(I_1/I_0)/dx\f$ 
*/
double h(double x){
	double tmp = (1.0 + ratio_bessel_20(x) - 2*SQ(ratio_bessel_10(x)) ); // besselI1_I0 replaced
	return(0.5*tmp);
}






/**
* @brief	Hessian matrix(w.r.t. W)
* @param 	W is a nx3 matrix.
* @param 	Psi_inv
* @param 	beta 
* @param 	TE 
* @param 	TR
* @param 	sigma
* @param 	r
* @param 	MRF_obj
* @param 	with_MRF = 1
* @param 	verbose = 1
* @return	Hessian Matrix
* @todo		Parallize the big loop
*
* @details 	Hessian_Matrix would be a sparse 3n x 3n matrix
* The values are: \f$d^2 l* / dW_{i'k'} dW_{ik}\f$
* {i,k} are arranged in a fashion so that k remains like a subdivision under division of i
* i.e., 
* l-th column/row representste k, i corresponding to : 
	k = l % 3; 
	i = (int) l/3;
* and accordingly,
	l = 3 * i + k;
* 
* Number of non-zero elements per column: <= 7*3
*/
SpMat Hessian_mat(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                  const Vector_eig &TE, const Vector_eig &TR, 
                  const Vector_eig &sigma, const Matrix_eig_row &r, 
                  MRF_param &MRF_obj, 
                  const Eigen::Matrix<char, Eigen::Dynamic, 1> &black_list, 
                  int with_MRF = 1, int verbose = 1){

	
	auto time_1_hess = std::chrono::high_resolution_clock::now();
	Debug1("Hessian calculation started");
	
	
	Matrix_eig_row v = v_mat(W, TE, TR);
	int n = W.rows();
	int m = v.cols();
	double temp = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;
	SpMat Gamma_inv;
	SpMat W_hess(3*n, 3*n);
	
	
	
	if(with_MRF){
		Gamma_inv = MRF_obj.Lambda(beta);
		W_hess.reserve( Eigen::VectorXi::Constant(3*n, 7*3) );
		Debug1("Hessian matrix allocated");
		
		// First, the Kroneker product term:
		SpMat Psi_inv_sp = Psi_inv.sparseView();
		W_hess = -Kron_Sparse_eig(Gamma_inv, Psi_inv_sp);
		Debug1("MRF part done of Hessian!");
	} else {
		W_hess.reserve( Eigen::VectorXi::Constant(3*n, 3) );
		Debug1("Hessian matrix allocated");
	}
	
	
	
	
	
	
	// Block-Diagonal parts //
	int i = 0, k = 0, k1 = 0, j = 0;
	Vector_eig temp_vec(3), temp_vec_1(3), temp_vec_2(3);;
	
	
	
//	Eigen::Matrix<char, Eigen::Dynamic, 1> black_list = Eigen::Matrix<char, Eigen::Dynamic, 1>::Ones(n);
//	
//	for(int i = 0; i < n; ++i){
//		for(int j = 0; j < m; ++j){
//			if(r(i, j) > 50){
//				black_list(i) = 0;
//				break;
//			}
//		}
//	}
//	if(verbose)
//		Debug0("Number of possible background voxels: " << (black_list.sum()));

	
	
	Matrix_eig A = Matrix_eig::Zero(3, 3);
	Eigen::LDLT<Matrix_eig> ldltOfA; // compute the Cholesky decomposition of A
	
	
	// Make it parallel
	for(i = 0; i < n; ++i){
		if(black_list(i) == 0){
		
			if(i % 10000 == 0){
				Rcpp::Rcout << "\n";
				Debug1("Hessian matrix, i: "<< i);
			}
			
			
			
			A.setZero(3, 3);
			for(k = 0; k < 3; ++k) {
				for(k1 = k; k1 < 3; ++k1) {
					
					temp = 0.;
					for(j = 0; j < m ; ++j) {
						
						tmp2 = r(i,j)/SQ(sigma(j));
						tmp3 = - v(i,j)/SQ(sigma(j)) + tmp2 * besselI1_I0(tmp2 * v(i,j));
						temp += tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1);
						
						tmp4 = (-1)/SQ(sigma(j)) + SQ(tmp2) * h(tmp2*v(i, j));
						
						temp += tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
										simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1);
						
						
						if(k == k1 && temp > 0.1){
							Debug1("i: " << i << ", j: " << j << ", k: " << k);
							Debug1("W.row(i): " << W.row(i) << "\t r.row(i): " << r.row(i) << "\tsigma(j): " << sigma(j));
							Debug1("tmp3: " << tmp3 << "\t tmp4: " << tmp4);
							Debug1("Added: 1st part: " << 
									tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1)
									<< ", 2nd part: " << 
									tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
										simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1));
							Debug1("final: " << (tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1) + 
										tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
											simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1)) << "\n");
							if(temp > 100){
								Debug0("Very Bad value in Hessian Matrix!!");
							}
						}
					}
					
					
					if(k == k1){
						A(k, k) = temp;
					} else {
						A(k, k1) = temp;
						A(k1, k) = temp;
					}
					
					// Chek the block is pd or not, then add
				}
			}
			
			
			ldltOfA.compute(A);
			if(ldltOfA.info() == Eigen::NumericalIssue){
				Rcpp::Rcout << "Not psd\n";
			} else {
				if(with_MRF){
					for(k = 0; k < 3; ++k) {
						for(k1 = k; k1 < 3; ++k1) {
							if(k == k1){
								W_hess.coeffRef(3 * i + k, 3 * i + k) += A(k, k);
							} else {
								W_hess.coeffRef(3 * i + k, 3 * i + k1) += A(k, k1);
								W_hess.coeffRef(3 * i + k1, 3 * i + k) += A(k, k1);
							}
						}
					}
				} else {
					for(k = 0; k < 3; ++k) {
						for(k1 = k; k1 < 3; ++k1) {
							if(k == k1){
								W_hess.insert(3 * i + k, 3 * i + k) = A(k, k);
							} else {
								W_hess.insert(3 * i + k, 3 * i + k1) = A(k, k1);
								W_hess.insert(3 * i + k1, 3 * i + k) = A(k, k1);
							}
						}
					}
				}
			}
		}
	}
	
	
	
	W_hess.makeCompressed();
	
	auto time_2_hess = std::chrono::high_resolution_clock::now();
	auto duration_hess = std::chrono::duration_cast<std::chrono::seconds>(time_2_hess - time_1_hess);
	Debug1("Time taken total loop: " << duration_hess.count() << " seconds\n");
	Debug0("Hessian calculated with MRF");
	
	return (-W_hess);	//3nx3n
}






/**
* @brief 	Calculates the \nu_ij for \nu_ij' \Sigma_ij \nu_ij (i.e., w.r.t. \nu_ij)
			where Sigma is an estimation of variance of (W_ik)_{i,k} 
 			i.e., it calculates \f$d\nu_{ij}/dW_{ik}\f$
* where j is fixed when we consider just one image and i corresponds to the i-th voxel
* 
* So, it would give a 3n x 1 vector: i.e., d\nu_ij/dW_{i1,k} (confusing notation)
* where the value is 0 when i != i1
* 
* So, we would get a sparse vector(/matrix) with 3 non-zero element of total size 3n x 1
* 
* 0:3 + 3*i - th elements would be non-zero and they are :
* d\nu_ij/dW_{i,0}, d\nu_ij/dW_{i,1}, d\nu_ij/dW_{i,2}
* 
* 
* grad is changed.
* 
* // https://stackoverflow.com/questions/47694725/using-inneriterator

  @param 	W
  @param 	Psi_inv
  @param 	beta
  @param 	TE
  @param 	TR
  @param 	sigma
  @param 	r
  @param 	i
  @param 	j
  @param[in, out] 	grad
  @return	Void
* 
*/
void v_grad(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
             const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, const Matrix_eig_row &r, 
             int i, int j, SpVec &grad){


	//SpMat grad(3*n_x*n_y*n_z, 1);
	//grad.reserve(VectorXi::Constant(1, 3));
	
	//for(int i1 = 3*i; i1 < 3 * i + 3; ++i1){
	//	grad.insert(i1, 1) = simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, i1 % 3);
	//}

	// SpVec grad(3*n_x*n_y*n_z);
	// Allocating this large vector might be cumbersome
	grad.setZero();
	
	
	for(int i1 = 3*i; i1 < 3 * i + 3; ++i1){
		// grad(i1) = simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, i1 % 3);
		grad.insert(i1) = simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, i1 % 3);
	}
	// return grad;
}











/**
* Hessian matrix with respect to new test images
* i.e., \nu_ij' \Sigma_ij^{-1} \nu_ij
* Iterative solution is used.

	@brief Variance corresponding to Test images using Information Matrix and Delta method
	@param	W			The W matrix
	@param 	Psi_inv
	@param 	beta 
	@param 	TE_train
	@param 	TR_train
	@param 	sigma_train
	@param 	train
	@param 	MRF_obj
	@param 	TE_test
	@param 	TR_test
	@param 	sigma_test
	@param 	test
	@param 	cg_maxiter
	@param 	cg_tol
	@param 	with_MRF = 1
	@return	Variance corresponding to the test images
*/
Matrix_eig Var_est_test_mat(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                            const Vector_eig &TE_train, const Vector_eig &TR_train,
                            const Vector_eig &sigma_train, const Matrix_eig_row &train, 
                            MRF_param &MRF_obj,
                            const Vector_eig &TE_test, const Vector_eig &TR_test, 
                            const Vector_eig &sigma_test, const Matrix_eig_row &test,
                            const Eigen::Matrix<char, Eigen::Dynamic, 1> &black_list,
                            int cg_maxiter = 1000, double cg_tol = 1e-6,
                            int with_MRF = 1){

	auto hess_1 = std::chrono::high_resolution_clock::now();
	int n = W.rows();
	Matrix_eig Var_est(n, TE_test.size());
	
	SpMat A = Hessian_mat(W, Psi_inv, beta, TE_train, TR_train, sigma_train, train, MRF_obj, black_list, with_MRF);
	assert(A.rows() == 3*n);
	// save_sparse(A, "Hessian_Matrix.csv", 1);
	
	
	
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
	Debug0("Diagonal");
//	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
//	Debug0("IncompleteCholesky");
	
	
	cg.setMaxIterations(cg_maxiter);
	cg.setTolerance(cg_tol);
	cg.compute(A);
	
	
	
	
	Vector_eig tmp_soln(3*n);
	SpVec b(3*n);
	// Vector_eig b = Vector_eig::Zero(3*n);
	
	
	// First i or j, which would be better? - check.
	Rcpp::Rcout << std::endl << std::flush;
	Debug1("Hessian loop starts!");
	Rcpp::Rcout << std::endl << std::flush;
	
	for(int j = 0; j < TE_test.size(); ++j){
		for(int i = 0; i < n; ++i){
		
			if( i % 100 == 0){
				Rcpp::Rcout << std::endl;
				Debug1("Info i: "<< i << ", j: " << j);
			}
			
			v_grad(W, Psi_inv, beta, TE_test, TR_test, sigma_test, test, i, j, b);
			
			tmp_soln = cg.solve(b);
			if( i % 1000 == 0)
				Rcpp::Rcout << "CG: #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
			Var_est(i, j) = b.dot(tmp_soln);			
		}
	}
	auto hess_2 = std::chrono::high_resolution_clock::now();
	auto hess_duration = std::chrono::duration_cast<std::chrono::seconds>(hess_2 - hess_1);
	Debug1("Time taken for Info matrix using Hessian: " << hess_duration.count() << " seconds\n");
	
	return Var_est;
}







/** @brief Variance corresponding to Contrast using Information Matrix and Delta method
	@param	W			The W matrix
	@param 	Psi_inv
	@param 	beta 
	@param 	TE_train
	@param 	TR_train
	@param 	sigma_train
	@param 	train
	@param 	MRF_obj
	@param 	TE_test
	@param 	TR_test
	@param 	sigma_test
	@param 	test
	@param 	contrast
	@param 	cg_maxiter
	@param 	cg_tol
	@param 	with_MRF = 1
	@return	Variance corresponding to the contrast vector
*/
Vector_eig Var_est_test_mat_contrast(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                                     const Vector_eig &TE_train, const Vector_eig &TR_train,
                                     const Vector_eig &sigma_train, const Matrix_eig_row &train, 
                                     MRF_param &MRF_obj,
                                     const Vector_eig &TE_test, const Vector_eig &TR_test, 
                                     const Vector_eig &sigma_test, const Matrix_eig_row &test,
                                     const SpVec &contrast, 
                                     const Eigen::Matrix<char, Eigen::Dynamic, 1> &black_list, 
                                     int cg_maxiter = 1000, double cg_tol = 1e-6,
                                     std::string preconditioner = "diagonal", 
                                     //auto preconditioner_2 = Eigen::DiagonalPreconditioner<double>, 
                                     int with_MRF = 1){

	auto hess_1 = std::chrono::high_resolution_clock::now();
	int n = W.rows();
	Matrix_eig Var_est(contrast.nonZeros(), TE_test.size());
	Vector_eig Var_est_vec(TE_test.size());
	
	
	SpMat A = Hessian_mat(W, Psi_inv, beta, TE_train, TR_train, sigma_train, train, MRF_obj, black_list, with_MRF);
	assert(A.rows() == 3*n);
	// save_sparse(A, "Hessian_Matrix.csv", 1);
	
	/*
	if(preconditioner == "diagonal"){
		Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
		Debug0("Diagonal");
	} else if(preconditioner == "incompletecholesky"){
		Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
		Debug0("IncompleteCholesky");
	}
	*/
	
	/*
	(preconditioner == "diagonal") ? Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
	: Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
	*/
	
	/*
	Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, 
		(preconditioner == "diagonal" ? Eigen::DiagonalPreconditioner<double> : Eigen::IncompleteCholesky<double>)> cg;
	*/
	
	/*
	Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, preconditioner_2> cg;
	*/
	
	Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
	Debug0("Diagonal");
	//Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
	//Debug0("IncompleteCholesky");
	
	
	cg.setMaxIterations(cg_maxiter);
	cg.setTolerance(cg_tol);
	
	cg.compute(A);
	Debug1("CG compute part done");
	
	
	
	
	Vector_eig tmp_soln(3*n);
	SpVec b(3*n);
	
	
	Vector_eig x = Vector_eig::Zero(3*n);
	for(int j = 0; j < TE_test.size(); ++j){
	
		x.setZero(3*n);
		// setzero
		for(SpVec::InnerIterator i_(contrast); i_; ++i_){
			v_grad(W, Psi_inv, beta, TE_test, TR_test, sigma_test, test, i_.index(), j, b);
			x += i_.value() * b;
		}
		
		auto hess_1_tmp = std::chrono::high_resolution_clock::now();		
		tmp_soln = cg.solve(x);
		auto hess_2_tmp = std::chrono::high_resolution_clock::now();
		auto hess_duration_tmp = std::chrono::duration_cast<std::chrono::seconds>(hess_2_tmp - hess_1_tmp);
		Debug1("Time taken for Info matrix using Hessian: " << hess_duration_tmp.count() << " seconds\n");
		Rcpp::Rcout << "CG: #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
		Var_est_vec(j) = x.dot(tmp_soln);
	}
	
	
	auto hess_2 = std::chrono::high_resolution_clock::now();
	auto hess_duration = std::chrono::duration_cast<std::chrono::seconds>(hess_2 - hess_1);
	Debug1("Time taken for Info matrix using Hessian: " << hess_duration.count() << " seconds\n");
	
	return Var_est_vec;
}






/** @brief Variance corresponding to New images using Parametric Bootstrap
	@param	W			The W matrix
	@param 	Psi_inv
	@param 	beta 
	@param 	TE_train
	@param 	TR_train
	@param 	sigma_train
	@param 	train
	@param 	r_scale
	@param 	TE_scale
	@param 	TR_scale
	@param 	MRF_obj
	@param 	TE_test
	@param 	TR_test
	@param 	sigma_test
	@param 	test
	@param 	B			Bootstrap number of replication
	@param 	EM_iter		Number of EM iterations
	@param 	abs_diff	absolute diff (parameters for the AECM)
	@param 	rel_diff	relative diff (parameters for the AECM)
	@param 	with_MRF = 1
	@return	Variance corresponding to the test images
*/
Matrix_eig para_boot_test_mat(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                              const Vector_eig &TE_train, const Vector_eig &TR_train, 
                              const Vector_eig &sigma_train, const Matrix_eig_row &r,
                              double r_scale, double TE_scale, double TR_scale, MRF_param &MRF_obj, 
                              const Vector_eig &TE_test, const Vector_eig &TR_test, 
                              const Vector_eig &sigma_test, const Matrix_eig_row &test,
                              const Eigen::Matrix<char, Eigen::Dynamic, 1> &black_list,
                              int B = 15, int EM_iter = 10, double abs_diff = 1.0e-4, double rel_diff = 1e-3,
                              int with_MRF = 1){

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
		AECM_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, generated_r, 
							r_scale, TE_scale, TR_scale, MRF_obj, black_list, EM_iter, with_MRF, abs_diff, rel_diff, 0);
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






/** @brief Variance corresponding to Contrast using Parametric Bootstrap
	@param	W			The W matrix
	@param 	Psi_inv
	@param 	beta 
	@param 	TE_train
	@param 	TR_train
	@param 	sigma_train
	@param 	train
	@param 	r_scale
	@param 	TE_scale
	@param 	TR_scale
	@param 	MRF_obj
	@param 	TE_test
	@param 	TR_test
	@param 	sigma_test
	@param 	test
	@param 	contrast
	@param 	B			Bootstrap number of replication
	@param 	EM_iter		Number of EM iterations
	@param 	abs_diff	absolute diff (parameters for the AECM)
	@param 	rel_diff	relative diff (parameters for the AECM)
	@param 	with_MRF = 1
	@return	Variance corresponding to the contrast vector
*/
Vector_eig para_boot_test_mat_contrast(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                                       const Vector_eig &TE_train, const Vector_eig &TR_train, 
                                       const Vector_eig &sigma_train, const Matrix_eig_row &r,
                                       double r_scale, double TE_scale, double TR_scale, MRF_param &MRF_obj, 
                                       const Vector_eig &TE_test, const Vector_eig &TR_test, 
                                       const Vector_eig &sigma_test, const Matrix_eig_row &test,
                                       const SpVec &contrast, 
                                       const Eigen::Matrix<char, Eigen::Dynamic, 1> &black_list,
                                       int B = 15, int EM_iter = 10, double abs_diff = 1e-1, double rel_diff = 1e-4, 
                                       int with_MRF = 1){

	auto boot_1 = std::chrono::high_resolution_clock::now();
	Debug1("Parametric Bootstrap starts!!");
	
	int n = W.rows();
	int m_test = TE_test.size();
	int m_train = TE_train.size();
	
	Vector_eig Var_est(m_test);
	Matrix_eig_row W_init = W;
	Matrix3d_eig Psi_inv_init = Psi_inv;
	Vector_eig beta_init = beta;
	
	// All estimated parametrs:
	Matrix_eig generated_r(n, m_train);		// train columns only
	
	// added matrices and vectors: 
	Matrix_eig tmp_mat = Matrix_eig::Zero(n, m_test);
	Vector_eig tmp_vec = Vector_eig::Zero(m_test);
	Vector_eig sum_mat = Vector_eig::Zero(m_test);
	Vector_eig sum_sq_mat = Vector_eig::Zero(m_test);
	
	for(int b = 0; b < B; ++b){
		if( b % 50 == 0){
			Debug0("bootstrap sample " << b);
		}

		
		//Generate an image matrix:
		generated_r = Gen_r(W, TE_train, TR_train, sigma_train);
		// W_init.noalias() = W;		// Not needed? - numerical stabilty?
		AECM_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, generated_r, 
							r_scale, TE_scale, TR_scale, MRF_obj, black_list, EM_iter, with_MRF, abs_diff, rel_diff, 0);
		tmp_mat = v_mat(W_init, TE_test, TR_test);
		// This is the nu_hat matrix for b-th replication - would be of  size n x m_test
		
		tmp_vec = tmp_mat.transpose() * contrast;
		sum_mat += tmp_vec;
		sum_sq_mat += tmp_vec.array().square().matrix();
	}
	
	sum_mat /= B;
	sum_sq_mat /= B;
	sum_sq_mat -= sum_mat.array().square().matrix();
	
	auto boot_2 = std::chrono::high_resolution_clock::now();
	auto boot_duration = std::chrono::duration_cast<std::chrono::seconds>(boot_2 - boot_1);
	Debug1("Time taken for Variance using Parametric Bootstrap: " << boot_duration.count() << " seconds\n");
	
	return (sum_sq_mat);
}





#endif /* !_VAR_HEADER_ */
