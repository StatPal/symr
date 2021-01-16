/* _VAR_HEADER_ */
#ifndef _VAR_HEADER_
#define _VAR_HEADER_




#include "functions_gen.hpp"
#include "functions_LS_and_init_value.hpp"
#include "functions_AECM.hpp"

#include <ctime>
#include <iomanip>









/***************************************************
**************** Information Matrix ****************
****************************************************/


/**
* Derivative of I_1/I_0
*/
double h(double x){
	double tmp = (1.0 + ratio_bessel_20(x) - 2*SQ(ratio_bessel_10(x)) ); // besselI1_I0 replaced
	return(0.5*tmp);
}






/**
* Hessian matrix(w.r.t. W)
* W is a nx3 matrix.
* Hessian_Matrix would be a sparse 3n x 3n matrix
* The values are: d^2 l* / dW_{i'k'} dW_{ik}
* {i,k} are arranged in a fashion so that k remains like a subdivision under division of i
* i.e., 
* l-th column/row representste k, i corresponding to : 
	k = l%3; 
	i = (int) l/3;
* and accordingly,
	l = 3 * i + k;
* 
* Number of non-zero elements per column: <= 7*3
*/
SpMat Hessian_mat(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                  const Vector_eig &TE, const Vector_eig &TR, 
                  const Vector_eig &sigma, const Matrix_eig_row &r, 
                  int n_x, int n_y, int n_z, MRF_param &MRF_obj, int with_MRF = 1, int verbose = 1){

	
	auto time_1_hess = std::chrono::high_resolution_clock::now();
	Debug1("Hessian calculation started");
	
	
	Matrix_eig_row v = v_mat(W, TE, TR);
	int n = n_x * n_y * n_z;
	int m = v.cols();
	double temp = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp31 = 0.0, tmp4 = 0.0;
	SpMat Gamma_inv;
	if(with_MRF){
		Gamma_inv = MRF_obj.Lambda(beta);
	}
	SpMat W_hess(3*n, 3*n);
	if(with_MRF){
		W_hess.reserve( Eigen::VectorXi::Constant(3*n, 7*3) );
	} else {
		W_hess.reserve( Eigen::VectorXi::Constant(3*n, 3) );
	}
	Debug1("Hessian matrix allocated");
	
	
	
	// First, the Kroneker prod term:
	if(with_MRF){
		SpMat Psi_inv_sp = Psi_inv.sparseView();
		W_hess = -Kron_Sparse_eig(Gamma_inv, Psi_inv_sp);
		Debug0("MRF part done of Hessian!");
	}
	
	
	
	
	// Diagonal parts //
	int i = 0, i1 = 0, k = 0, k1 = 0, j = 0;
	Vector_eig temp_vec(3), temp_vec_1(3), temp_vec_2(3);;
	
	
	for(i = 0; i < n; ++i) {
	
		if(i % 10000 == 0){
			std::cout << "\n";
			Debug1("Hess matrix i: "<< i << ", j: " << j);
		}
		
		
		//temp_vec = W.row(i);
		for(k = 0; k < 3; ++k) {
			//for(k1 = 0; k1 < 3; ++k1) {
			for(k1 = k; k1 < 3; ++k1){
				
				temp = 0.;
				for(j = 0; j < m ; ++j) {
					
					tmp2 = r(i,j)/SQ(sigma(j));
					//if(i == 0){
						//Debug0(tmp2);
					//}
					tmp3 = - v(i,j)/SQ(sigma(j)) + tmp2 * besselI1_I0(tmp2 * v(i,j));
					//if(i == 0){
						//Debug0(tmp3);
					//}
					temp += tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1);
					//if(i == 0){
						//Debug0(simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1));  // problem
						//Debug0(tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1));
					//}
					
					
					//tmp2 *= v(i,j);
					//tmp3 = (1 + ratio_bessel_20(tmp2) - 2*SQ(besselI1_I0(tmp2)) );
					//tmp4 = -1/SQ(sigma(j)) +  0.5*SQ(r(i,j)/SQ(sigma(j)))*tmp3;
					// This is also valid
					
					
					
					
					// tmp4 = (-1)/SQ(sigma(j)) + SQ(tmp2) * h(tmp2*v(i, j)/SQ(sigma(j)));  // BUG
					tmp4 = (-1)/SQ(sigma(j)) + SQ(tmp2) * h(tmp2*v(i, j));
					// This is also valid
					
					//if(i == 0){
					//	Debug0(tmp2);
					//	Debug0(v(i,j));
					//	Debug0(SQ(sigma(j)));
					//	Debug0(tmp2*v(i, j)/SQ(sigma(j)));
					//	Debug0(h(tmp2*v(i, j)/SQ(sigma(j))));
					//	Debug0(tmp4);
					//}
					temp += tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
									simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1);
					//if(i == 0){
					//	Debug0(simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k));
					//	Debug0(simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
					//				simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1));
					//	Debug0(tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
					//				simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1));
					//	Debug0("Added:");
					//	Debug0(tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1));
					//	Debug0(tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
					//				simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1));
					//	Debug0((tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1) + 
					//				tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
					//					simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1)));
					//	Debug0("\n")
					//}
					//if(i == 0){
					//	Debug0(temp);
					//	Debug0("\n\n");
					//}
					
					if(k == k1 && temp > 0.1){
						Debug0("i: " << i << ", j: " << j << ", k: " << k);
						Debug0("W.row(i): " << W.row(i) << "\t r.row(i): " << r.row(i) << "\tsigma(j): " << sigma(j));
						Debug0("tmp3: " << tmp3 << "\t tmp4: " << tmp4);
						Debug0("Added: 1st part: " << tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1)
								<< ", 2nd part: " << tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
									simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1));
						Debug0("final: " << (tmp3 * simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W.row(i), TE, TR, j, k, k1) + 
									tmp4 * simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k) * 
										simple_dee_v_ij_dee_W_ik(W.row(i), TE, TR, j, k1)) << "\n");
						if(temp > 100){
							Debug0("very high!!");
							// exit(EXIT_FAILURE);
							// change
						}
					}
					
				}
				// W_hess.insert(i+k*n, i+k1*n) = temp;		// old way - not very visually pleasing I guess.
				
				if(with_MRF){
					if(k == k1){
						W_hess.coeffRef(3 * i + k, 3 * i + k) += temp;
						// BUG? check negativity and positivity
						// minus added for Gamma * Psi						
					} else {
						W_hess.coeffRef(3 * i + k, 3 * i + k1) += temp;
						W_hess.coeffRef(3 * i + k1, 3 * i + k) += temp;
					}
				} else {
					if(k == k1){
						W_hess.insert(3 * i + k, 3 * i + k) = temp;
					} else {
						W_hess.insert(3 * i + k, 3 * i + k1) = temp;
						W_hess.insert(3 * i + k1, 3 * i + k) = temp;
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






// https://stackoverflow.com/questions/47694725/using-inneriterator
/**
* Calculates the \nu_ij for \nu_ij' \Sigma_ij \nu_ij (i.e., w.r.t. \nu_ij)
* where Sigma is an estimation of variance of (W_ik)_{i,k} 
* 
* i.e., it calculates d\nu_ij/dW_ik
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
*/
void v_grad(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
             const Vector_eig &TE, const Vector_eig &TR, const Vector_eig &sigma, const Matrix_eig_row &r, 
             int n_x, int n_y, int n_z, int i, int j, SpVec &grad){


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
	save_sparse(A, "Hessian_Matrix.csv", 1);
	
	// Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
	// Debug0("Diagonal");
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
	Debug0("IncompleteCholesky");
	
	
	cg.setMaxIterations(80);
//	cg.setTolerance(1e-10);
	cg.compute(A);
	
	
	
	
	Vector_eig tmp_soln(3*n);
	SpVec b(3*n);
	// Vector_eig b = Vector_eig::Zero(3*n);
	
	
	// First i or j, which would be better? - check.
	std::cout << std::endl << std::flush;
	Debug1("Hessian loop starts!");
	std::cout << std::endl << std::flush;
	for(int j = 0; j < TE_test.size(); ++j){
		for(int i = 0; i < n; ++i){
		
			if( i % 100 == 0){
				std::cout << std::endl;
				Debug1("Info i: "<< i << ", j: " << j);
			}
			
			// b = v_grad(W, Psi_inv, beta, TE_test, TR_test, sigma_test, test, n_x, n_y, n_z, i, j);
			v_grad(W, Psi_inv, beta, TE_test, TR_test, sigma_test, test, n_x, n_y, n_z, i, j, b);
			
			assert(b.size() == A.cols());
			
			tmp_soln = cg.solve(b);
			if( i % 100 == 0)
				std::cout << "CG: #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
			Var_est(i, j) = b.dot(tmp_soln);
			
			// b(3*i) = 0.0; b(3*i+1) = 0.0; b(3*i+2) = 0.0;
			
		}
	}
	auto hess_2 = std::chrono::high_resolution_clock::now();
	auto hess_duration = std::chrono::duration_cast<std::chrono::seconds>(hess_2 - hess_1);
	Debug1("Time taken for Info matrix using Hessian: " << hess_duration.count() << " seconds\n");
	
	return Var_est;
}







/**
* Variance corresponding to Contrast using Information Matrix and Delta method
*/
Vector_eig Var_est_test_mat_contrast(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
                                     const Vector_eig &TE_train, const Vector_eig &TR_train,
                                     const Vector_eig &sigma_train, const Matrix_eig_row &train, 
                                     int n_x, int n_y, int n_z, MRF_param &MRF_obj,
                                     const Vector_eig &TE_test, const Vector_eig &TR_test, 
                                     const Vector_eig &sigma_test, const Matrix_eig_row &test,
                                     const SpVec &contrast){

	auto hess_1 = std::chrono::high_resolution_clock::now();
	int n = W.rows();
	Matrix_eig Var_est(contrast.nonZeros(), TE_test.size());
	
	SpMat A = Hessian_mat(W, Psi_inv, beta, TE_train, TR_train, sigma_train, train, n_x, n_y, n_z, MRF_obj, 1);	
	//SpMat A = Hessian_mat(W, Psi_inv, beta, TE_train, TR_train, sigma_train, train, n_x, n_y, n_z, MRF_obj, 0);
	assert(A.rows() == 3*n);
	save_sparse(A, "Hessian_Matrix.csv", 1);
	
	// Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
	// Debug0("Diagonal");
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
	Debug0("IncompleteCholesky");
	cg.compute(A);
	
	
	
	
	
	Vector_eig tmp_soln(3*n);
	SpVec b(3*n);
	// Vector_eig b = Vector_eig::Zero(3*n);
	
	
	// First i or j, which would be better? - check.
	std::cout << std::endl << std::flush;
	Debug1("Hessian loop starts!");
	std::cout << std::endl << std::flush;
	for(int j = 0; j < TE_test.size(); ++j){
	
		for(SpVec::InnerIterator i_(contrast); i_; ++i_){
		
			Debug1("Info i: "<< i_.index() << ", j: " << j);
			
			// b = v_grad(W, Psi_inv, beta, TE_test, TR_test, sigma_test, test, n_x, n_y, n_z, i_.index(), j);
			v_grad(W, Psi_inv, beta, TE_test, TR_test, sigma_test, test, n_x, n_y, n_z, i_.index(), j, b);
			
			assert(b.size() == A.cols());
			
			tmp_soln = cg.solve(b);
			std::cout << "CG: #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
			Var_est(i_.index(), j) = b.dot(tmp_soln);
			
			// b(3*i_.index()) = 0.0; b(3*i_.index()+1) = 0.0; b(3*i_.index()+2) = 0.0;
			
			Var_est(i_.index(), j) *= SQ(i_.value()); 
			
		}
	}
	
	
	
	auto hess_2 = std::chrono::high_resolution_clock::now();
	auto hess_duration = std::chrono::duration_cast<std::chrono::seconds>(hess_2 - hess_1);
	Debug1("Time taken for Info matrix using Hessian: " << hess_duration.count() << " seconds\n");
	
	return Var_est.array().colwise().sum();
}






/**
* Parametric Bootstrap
* for test set of images
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
		AECM_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, generated_r, 
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
