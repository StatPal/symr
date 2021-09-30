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
#ifndef _VAR_HEADER_3D_
#define _VAR_HEADER_3D_




#include "../functions_gen_VAR.hpp"
#include "functions_AECM.hpp"

#include <ctime>
#include <iomanip>




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
Matrix_eig para_boot_test_mat_3D(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
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
Vector_eig para_boot_test_mat_contrast_3D(const Matrix_eig_row &W, const Matrix3d_eig &Psi_inv, const Vector_eig &beta, 
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





#endif /* !_VAR_HEADER_3D_ */
