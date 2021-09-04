//#include <RcppEigen.h>
//#include <RcppGSL.h>

// c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) AECM_wrapper.cpp -o sympy$(python3-config --extension-suffix) -I /usr/include/eigen3 -lgsl -lgslcblas -lm -fopenmp



#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

using namespace pybind11::literals;



#include "include/functions_gen.hpp"
#include "include/functions_LS_and_init_value.hpp"

#include "include/2D/functions_AECM.hpp"
#include "include/2D/functions_OSL.hpp"



typedef Eigen::MappedSparseMatrix< double > mappedSparseMatrix ;
//typedef Eigen::Map<Eigen::MatrixXd> mappedMatrix;
//typedef Eigen::Map<Eigen::VectorXd> mappedVector;


// add mean_rice, J_n matrix and eigenvalues
// Bloch vector, v_mat, generation of r matrix
// Derivatives and double derivatves
// Hessian etc etc etc


double mean_rice_R(double nu, double sigma){
	return mean_rice(nu, sigma);
}


Eigen::SparseMatrix<double> J(int n){
	return J_n(n);
}


Eigen::VectorXd eigenvals_J(int n) {
	return eigenvals_J_n(n);
}


Eigen::VectorXd Bloch_eqn_R(const Eigen::VectorXd W_row, 
							const Eigen::VectorXd TE, const Eigen::VectorXd TR){
	Eigen::VectorXd tmp = Eigen::Map<Eigen::VectorXd>::Zero(TE.size());
	Bloch_vec(W_row, TE, TR, tmp);
	return tmp;
}


Eigen::MatrixXd v_mat_R(Eigen::MatrixXd &W, 
						Eigen::VectorXd &TE, Eigen::VectorXd &TR){
	return v_mat(W, TE, TR);
}


Eigen::MatrixXd Generate_r(const Eigen::MatrixXd &W, 
							const Eigen::VectorXd &TE, const Eigen::VectorXd &TR, 
					     	const Eigen::VectorXd &sigma){
	return Gen_r(W, TE, TR, sigma);
}



double dee_v_ij_dee_W_ik(const Eigen::VectorXd &W_row, 
						const Eigen::VectorXd &TE, const Eigen::VectorXd &TR, 
						int j, int k){
	return simple_dee_v_ij_dee_W_ik(W_row, TE, TR, j-1, k-1);
}


double dee_2_v_ij_dee_W_ik_dee_W_ik1(const Eigen::VectorXd &W_row, 
									const Eigen::VectorXd &TE, const Eigen::VectorXd &TR, 
                                    int j, int k, int k1){
	return simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W_row, TE, TR, j-1, k-1, k1-1);
}





Eigen::MatrixXd Init_val_least_sq_R(const Eigen::MatrixXd &train, 
                        const Eigen::VectorXd &TE_train, const Eigen::VectorXd &TR_train, 
                        Eigen::VectorXd our_dim_1, 
                        double train_scale, double TE_scale, double TR_scale, int maxiter_LS, 
                        double W_1_init, double W_2_init){ //, 
                        // int do_least_sq = 1, char will_write = 0){
    
    
    //Primary Initial value for test//
	int n = train.rows();
	Matrix_eig_row W = Matrix_eig_row::Ones(train.rows(), 3);
	
	
	// Ad hoc initial values:
	W.col(0) = train.rowwise().mean().transpose();
	W.col(1) *= W_1_init;
	W.col(2) *= W_2_init;
	for(int i = 0; i < n; ++i){
		if(W(i, 0)>450.0){
			W(i, 0) = 425.0;
		}
		if(W(i, 0) < 0.0001){		// Added later to recover from sub lb[0] case possible due to scaling sown
			W(i, 0) = 0.0001;
		}
		
	}
	Debug1("1st level preprocess of initial value done!\n----------------\n----------------\n");
		
	least_sq_solve(W, TE_train, TR_train, train, train_scale, TE_scale, TR_scale, maxiter_LS);
	
	return(W);
}






#include <tuple>




std::tuple<Eigen::MatrixXd, Eigen::Matrix3d, Eigen::VectorXd> AECM_R(Eigen::MatrixXd W, Eigen::VectorXd our_dim_1,
            const Eigen::VectorXd &TE_train, const Eigen::VectorXd &TR_train, 
            const Eigen::VectorXd &sigma_train, const Eigen::MatrixXd &train, 
            double train_scale, double TE_scale, double TR_scale, 
            const Eigen::VectorXd &black_list, 
	        int maxiter = 50, int penalized = 1, 
            double abs_diff = 1e-1, double rel_diff = 1e-5, int verbose = 0, int verbose2 = 0){


	// Matrix_eig_row W_init = Matrix_eig_row::Zero(train.rows(), 3);	// Change this format and take initial value as input
	Matrix_eig_row W_init = W;
	
	short our_dim_train[8];
	for(int i = 0; i < 8; ++i){
		our_dim_train[i] = our_dim_1[i];
	}

	Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
	Vector_eig beta_init = 1.0*Vector_eig::Ones(3);
	
	MRF_param MRF_obj_1(our_dim_train[1], our_dim_train[2], our_dim_train[3]);
	
	Eigen::Matrix<char, Eigen::Dynamic, 1> black_list_2 = black_list.cast<char>();
	
	
	AECM_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
	          train_scale, TE_scale, TR_scale, MRF_obj_1, black_list_2, 
	          maxiter, penalized, abs_diff, rel_diff, verbose, verbose2);

	return std::make_tuple(W_init, Psi_inv_init, beta_init);
}










std::tuple<Eigen::MatrixXd, Eigen::Matrix3d, Eigen::VectorXd> OSL_R(Eigen::MatrixXd W, Eigen::VectorXd our_dim_1,
            const Eigen::VectorXd &TE_train, const Eigen::VectorXd &TR_train, 
            const Eigen::VectorXd &sigma_train, const Eigen::MatrixXd &train, 
            double train_scale, double TE_scale, double TR_scale, 
            const Eigen::VectorXd &black_list, 
	        int maxiter = 50, int penalized = 1, 
            double abs_diff = 1e-1, double rel_diff = 1e-5, int verbose = 0, int verbose2 = 0){


	// Matrix_eig_row W_init = Matrix_eig_row::Zero(train.rows(), 3);	// Change this format and take initial value as input
	Matrix_eig_row W_init = W;
	
	short our_dim_train[8];
	for(int i = 0; i < 8; ++i){
		our_dim_train[i] = our_dim_1[i];
	}

	Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
	Vector_eig beta_init = 1.0*Vector_eig::Ones(3);
	
	MRF_param MRF_obj_1(our_dim_train[1], our_dim_train[2], our_dim_train[3]);
	
	Eigen::Matrix<char, Eigen::Dynamic, 1> black_list_2 = black_list.cast<char>();
	
	
	OSL_optim(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
	          train_scale, TE_scale, TR_scale, MRF_obj_1, black_list_2, 
	          maxiter, penalized, abs_diff, rel_diff, verbose, verbose2);

	
	
	return std::make_tuple(W_init, Psi_inv_init, beta_init);
}






Eigen::VectorXd Performance_test_R(const Eigen::MatrixXd &W, const Eigen::MatrixXd &test, 
							const Eigen::VectorXd &TE_test, const Eigen::VectorXd &TR_test, 
							const Eigen::VectorXd &sigma_test,const Eigen::VectorXd &black_list, 
							int v_type = 1, int measure_type = 1, int scale = 1, int verbose = 0){


	Eigen::Matrix<char, Eigen::Dynamic, 1> black_list_2 = black_list.cast<char>();
	
	return(Performance_test(W, test, TE_test, TR_test, sigma_test, black_list_2,
					 v_type, measure_type, scale, verbose));
}












// 3D

#include "include/3D/functions_AECM.hpp"
#include "include/3D/functions_OSL.hpp"




std::tuple<Eigen::MatrixXd, Eigen::Matrix3d, Eigen::VectorXd> AECM_R_3D(Eigen::MatrixXd W, Eigen::VectorXd our_dim_1,
            const Eigen::VectorXd &TE_train, const Eigen::VectorXd &TR_train, 
            const Eigen::VectorXd &sigma_train, const Eigen::MatrixXd &train, 
            double train_scale, double TE_scale, double TR_scale, 
            const Eigen::VectorXd &black_list, 
	        int maxiter = 50, int penalized = 1, 
            double abs_diff = 1e-1, double rel_diff = 1e-5, int verbose = 0, int verbose2 = 0){


	// Matrix_eig_row W_init = Matrix_eig_row::Zero(train.rows(), 3);	// Change this format and take initial value as input
	Matrix_eig_row W_init = W;
	
	short our_dim_train[8];
	for(int i = 0; i < 8; ++i){
		our_dim_train[i] = our_dim_1[i];
	}

	Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
	Vector_eig beta_init = 1.0*Vector_eig::Ones(3);
	
	MRF_param MRF_obj_1(our_dim_train[1], our_dim_train[2], our_dim_train[3]);
	
	Eigen::Matrix<char, Eigen::Dynamic, 1> black_list_2 = black_list.cast<char>();
	
	
	AECM_optim_3D(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
	          train_scale, TE_scale, TR_scale, MRF_obj_1, black_list_2, 
	          maxiter, penalized, abs_diff, rel_diff, verbose, verbose2);

	return std::make_tuple(W_init, Psi_inv_init, beta_init);
}







std::tuple<Eigen::MatrixXd, Eigen::Matrix3d, Eigen::VectorXd> OSL_R_3D(Eigen::MatrixXd W, Eigen::VectorXd our_dim_1,
            const Eigen::VectorXd &TE_train, const Eigen::VectorXd &TR_train, 
            const Eigen::VectorXd &sigma_train, const Eigen::MatrixXd &train, 
            double train_scale, double TE_scale, double TR_scale, 
            const Eigen::VectorXd &black_list, 
	        int maxiter = 50, int penalized = 1, 
            double abs_diff = 1e-1, double rel_diff = 1e-5, int verbose = 0, int verbose2 = 0){


	// Matrix_eig_row W_init = Matrix_eig_row::Zero(train.rows(), 3);	// Change this format and take initial value as input
	Matrix_eig_row W_init = W;
	
	short our_dim_train[8];
	for(int i = 0; i < 8; ++i){
		our_dim_train[i] = our_dim_1[i];
	}

	Eigen::Matrix3d Psi_inv_init = Eigen::Matrix3d::Identity();
	Vector_eig beta_init = 1.0*Vector_eig::Ones(3);
	
	MRF_param MRF_obj_1(our_dim_train[1], our_dim_train[2], our_dim_train[3]);
	
	Eigen::Matrix<char, Eigen::Dynamic, 1> black_list_2 = black_list.cast<char>();
	
	
	OSL_optim_3D(W_init, Psi_inv_init, beta_init, TE_train, TR_train, sigma_train, train, 
	          train_scale, TE_scale, TR_scale, MRF_obj_1, black_list_2, 
	          maxiter, penalized, abs_diff, rel_diff, verbose, verbose2);

	return std::make_tuple(W_init, Psi_inv_init, beta_init);
}








// ----------------
// Python interface
// ----------------


PYBIND11_MODULE(sympy, m) {
    m.doc() = "pybind11 plugin corresponding to symr"; // optional module docstring

//    m.def("inv", &inv);
//    m.def("det", &det);
    
    m.def("mean_rice", 			&mean_rice_R);
    //m.def("J_n", 				&J_n);
    m.def("eigenvals_J_n", 		&eigenvals_J_n);
    m.def("Bloch_eqn", 			&Bloch_eqn_R);		// PROBLEM sarts here
    m.def("v_mat", 				&v_mat_R);
    m.def("Generate_r", 		&Generate_r);
    m.def("dee_v_ij_dee_W_ik", 		&dee_v_ij_dee_W_ik);
    m.def("dee_2_v_ij_dee_W_ik_dee_W_ik1", 		&dee_2_v_ij_dee_W_ik_dee_W_ik1);
    
    m.def("LS_est", 			&Init_val_least_sq_R);
    m.def("AECM_est",	 		&AECM_R);
    m.def("OSL_est", 			&OSL_R);
    m.def("AECM_est_3D",	 		&AECM_R_3D);
    m.def("OSL_est_3D", 			&OSL_R_3D);
    m.def("Performance_test", 	&Performance_test_R);
}


	


