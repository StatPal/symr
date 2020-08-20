/**
* Calculating estimates of \sigma_j's

g++ Var_3.cpp -o test_var -I /usr/include/eigen3 -O3

g++ Var_3.cpp -o test_var -I ../eigen-3.3.7 -O3



./test_var ../Read_Data/ZHRTS1.nii Dummy_sd.txt 0

./test_var ../Read_Data/new_phantom.nii Dummy_sd.txt 0


Init value my code
and EM translated from Dr. Maitra's code

*/



#include "scheme_new_numerical.hpp"
#include "Read_files_2.hpp"
#include "Init_value_6_numerical.hpp"

#include <fstream>
#include <random>

#define Inf 1e+140




// Initial value part:

int choose_rnd_prob(const Vector_eig &p){

	
	int n = p.size(), choose;
	// Normalize p first -- don't do this - makes the intervals very small and makes it hard to compare
	// Vector_eig p_norm = (p/p.sum()).eval();
	Vector_eig cum_sum_p = p; 			// p_norm;
	for(int i = 1; i < n; ++i){
		cum_sum_p(i) += cum_sum_p(i-1);
	}
	
	
	// Random no:
	//std::srand((unsigned int) time(0));	// https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	std::random_device rd{};
	std::mt19937 gen{rd()};
	std::uniform_real_distribution<> d{0.0, p.sum()}; // d{0.0,1.0};
	double random_val = d(rd);	//VectorXd::Random(1)
	
	if(0.0 <= random_val && random_val < cum_sum_p(0)){
		choose  = 0;
	}
	for(int i = 1; i < n; ++i){
		if(cum_sum_p(i - 1) <= random_val && random_val < cum_sum_p(i)){
			choose = i;
			break;
		}
	}
	
	return choose;
}


int closest_dist(double R_i, Vector_eig mu){

	Vector_eig diff = (Vector_eig::Constant(mu.size(), R_i) - mu).cwiseAbs();
	std::ptrdiff_t posn;
	int maxdiff = diff.minCoeff(&posn);
	int posn_1 = (int) posn;
	return (posn_1);
}



// Modified version:
Matrix_eig Est_var_init(const Vector_eig &R, int J, 
						double &sig_sq, Vector_eig &pi, Vector_eig &mu, int maxiter = 10){


	//Declare variables:
	int n = R.size();
	assert(pi.size() == J && mu.size() == J);
	Vector_eig W(n), X(n);
	sig_sq = 0.0;
	mu = Vector_eig::Constant(J, -100);
	// Numerical reasons - all (-100) - so that it's easier to calculate closest mu_j
	
	
	// Track the best:
	Vector_eig mu_best = mu, pi_best = pi;
	double sig_sq_best = sig_sq;
	
	
	
	// Step 1: - not repeated
	mu(J-1) = 0.0;		// or randomly chosen r
	W = Vector_eig::Constant(n, J-1);
	
	
	Vector_eig tmp = Vector_eig::Zero(n), diff(n);
	// Basically written as Step 4:
	int iter = 0;
	while(iter++ < maxiter){
		for(int j = 0; j < J-1; ++j){
		
			// std::cout << "\n";
			Debug2("j = " << j);
			
			// Step 2:
			for(int i = 0; i < n; ++i){
				X(i) = R(i) - mu(W(i));		// X can be -ve.
			}
			sig_sq = X.squaredNorm()/(2*n);
			Debug2("sig hat sq: " << sig_sq);
			
			
			// Step 3:
			int ind = choose_rnd_prob( 1 - exp( - X.array().pow(2.0) /(2*sig_sq)));
			// mu(j) = X(ind);
			// problem: X would be -ve also, creating -ve value of mu
			mu(j) = R(ind);					// changed 
			for(int i = 0; i < n; ++i){		// Getting the closest
				W(i) = closest_dist(R(i), mu);
			}
			
			Debug2("In the loop: mu size: " << mu.size() << " mu: " << mu.transpose() << "\n");
		}
		Debug2("Step 4 done!\n");
		
		// Step 5:
		Vector_eig mu_tmp = Vector_eig::Zero(J);
		pi = Vector_eig::Zero(J);
		for(int j = 0; j < J; ++j ){
			for(int i = 0; i < n; ++i){
				if(W(i) == j){
					pi(j) += (1.0);
					mu_tmp(j) += R(i);
				}
			}
			//Debug2("In the loop: mu_tmp: " << mu_tmp.size() << "; mu_tmp: " << mu_tmp.transpose());
			mu_tmp(j) /= pi(j);		// # of elements in that class
		}
		
		
		pi /= n;
		Debug2("Init pi: " << pi.transpose());
		
		mu_tmp(J-1) = 0.0;
		mu = mu_tmp;
		Debug2("Init mu: " << mu.transpose());
		
		for(int i = 0; i < n; ++i){
			X(i) = R(i) - mu(W(i));		// X can be -ve.
		}
		sig_sq = X.squaredNorm()/(2*n);
		Debug2("Init sig hat sq: " << sig_sq << "\nInitial value done!!\n\n" );
	}
	
	
	Matrix_eig W_mat = Matrix_eig::Zero(n, J);
	// Somehow, the mu - sig values are not doing good - so i am feeding hard - clustered W first
	for(int j = 0; j < J; ++j){
		for(int i = 0; i < n; ++i){
			W_mat(i, W(i)) = 1;
		}
	}
	
	Debug1("Init sig hat sq: " << sig_sq << "\nInit pi: " << pi.transpose() << ";\nInit mu: " << mu.transpose() << "\nInitial value done!!");
	
	
	
	//Checks:
	/*
	std::cout << "\n\n\n\n";
	Debug0("Checks!\n");
	Vector_eig tmp_mean = Vector_eig::Zero(J);
	Vector_eig tmp_sum = W_mat.colwise().sum();
	Vector_eig tmp_pi = W_mat.colwise().mean();
	double tmp_sig = 0.0;
	for(int j = 0; j < J; ++j){
		for(int i = 0; i < n; ++i){
			if(W_mat(i,j) == 1.0){
				tmp_mean(j) +=  R(i);
			}
		}
		tmp_mean(j) /= tmp_sum(j);
	}
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < J; ++j){
			if(W_mat(i,j) == 1.0){
				tmp_sig += SQ(R(i) - tmp_mean(j));
			}
		}
	}
	tmp_sig /= 2*n;
	Debug0("Check: tmp_pi: " << tmp_pi.transpose() << "; tmp_mean: " << tmp_mean.transpose() << 
			"; tmp_sig:" << tmp_sig);
	
	
	// Okay!
	*/
	
	return W_mat;			// Check this is good or not!
}















/**
log density of rice distn.
*/
double dlrice(double x, double mu, double sig_sq) {

	//double inf = std::numeric_limits<double>::infinity();
	if (x > 0.0) {
		double xv = x / sig_sq;
		if (mu > 0.0) {
			return (std::log(xv) - (SQ(x) + SQ(mu))/(2 * sig_sq) + logBesselI0(xv * mu));
		} else if (mu == 0.0){
			return (std::log(xv) - SQ(x)/(2 * sig_sq));
		} else {
			return -Inf;
		}
	} else {
		exit(EXIT_FAILURE);
	}
}










/*   This is the E-Step: assigns the responsibilities of each of the k groups 
     to each of the n observations. 
     The inputs are:
     n     = number of Ricean observations
     k     = number of components (all Rice; Rayleigh incorporated within, in 
             the last k'th component if present)
     X     = vector of Ricean observations
     pi    = vector of prior probabilities
     Mu    = vector of Ricean-distributed means (last zero, if Rayleigh 
             component is present)     
     Sigma = (common) noise parameter of Rice-distributed components
     Gamma = n x k matrix of posterior probability of observation X[i] 
             belonging to the kth group. 
     RGamm = n x k matrix of posterior expectation of the product of the
             indicator (that the ith observation is in kth group) and 
	     the cosine of angle given ith magnitude observation belongs to 
	     the kth group
*/

void rice_estep(int n, int k, const Vector_eig &X, Matrix_eig &Gamma, Matrix_eig &RGamm,
		Vector_eig &pi, Vector_eig &Mu, double Sigma){

	int i, l;
	double sum;
	
	Vector_eig temp = Vector_eig::Zero(k);
	Vector_eig minmax(2);
	
	for (i = 0; i < n; i++) {
		sum = 0.;
		for (l = 0; l < k; l++){
			temp[l] = std::log(pi[l]) + dlrice(X[i], Mu[l], Sigma);
			// I guess log(0) is creating problem in log pi case
			// Oh! mu and sigma were interchanged - why sigma is written before mu-confusing!!!!!!!!
			if(std::isnan(temp[l])){
				Debug1("l: "<< l << "dlrice(X[i], Mu[l], Sigma): " << dlrice(X[i], Mu[l], Sigma) << 
						"X[i], Sigma, Mu[l]:" << X[i] << ", " << Sigma << ", " <<  Mu[l]);
			}
		}
		
		
		//minmax = range(temp, k); /* get the minimum and maximum values*/
		minmax(0) = temp.array().minCoeff();
		minmax(1) = temp.array().maxCoeff();
		for (l = 0; l < k; l++) {
			temp[l] -= minmax[1]; /*reduce to get smaller values*/ 
			// Subtract maximum?? -Subrata  - check the fn - checked
			sum += std::exp(temp[l]);
		}
		
		for (l = 0; l < k; l++) {
			//Gamma[i][l] = std::exp(temp[l])/sum;
			//RGamm[i][l] = Gamma[i][l] * besselI1_I0( X[i] * Mu[l] / SQ(Sigma));
			Gamma(i, l)= std::exp(temp[l])/sum;
			RGamm(i, l) = Gamma(i, l) * besselI1_I0( X[i] * Mu[l] / SQ(Sigma));
		}
	}
}


void rice_mstep(const Vector_eig &X, int n, int k, Vector_eig &pi, Vector_eig &Mu, 
		double &Sigma, Matrix_eig &Gamma, Matrix_eig &RGamm){

	int i, ll, iRayleigh = 0;
	double sums = 0;
	
	Vector_eig sum = Vector_eig::Zero(k);

	if (Mu[k - 1] == 0) iRayleigh = 1; /*this means the Rayleigh component*/



	for (ll = 0; ll < k; ll++) {
		sum[ll]=0.;
		Mu[ll] = 0;
		for (i = 0; i < n; i++) {
			sum[ll] += Gamma(i, ll);
			Mu[ll] += X[i] * RGamm(i, ll);
		}
		Mu[ll] /= sum[ll];
		pi[ll] = sum[ll]/n;
	}
	
	
	if (iRayleigh) Mu[k - 1] = 0;
	
	for (i = 0; i < n; i ++) {
		sums += SQ(X[i]);
		for (ll = 0; ll < k; ll++) {
			sums -= 2 * Mu[ll] * X[i] * RGamm(i, ll);
			sums += Gamma(i, ll) * SQ(Mu[ll]);
		}
	}
	sums /= 2 * n;

	Sigma = sqrt(sums);
	
	Debug2("sums: " << sums << " Sigma:" << Sigma);
	
}



Eigen::VectorXi classify(int n, int k, Matrix_eig &Gamma);

double icl(int n, int k, const Vector_eig &X, Eigen::VectorXi &class1, 
			Vector_eig &Mu, Vector_eig &pi, double sigma);

double observedDataLogLikelihood(const Vector_eig &y, const int numPoints,
				 const Vector_eig &p, const Vector_eig &mu,
				 const int numClusters, const double sigma);



Eigen::VectorXi rice_emcluster(int n, int k, Vector_eig &pi, const Vector_eig &X, Vector_eig &Mu, 
								double &Sigma, int maxiter, double eps, double &llhdval,
								double &ICL, double &BIC){

	int iter = 0;
	Eigen::VectorXi class1(n);
	double oldllhd;
	
	Matrix_eig gamm = Matrix_eig::Zero(n, k);
	Matrix_eig rgamm = Matrix_eig::Zero(n, k);
	
	
	llhdval = observedDataLogLikelihood(X, n, pi, Mu, k, Sigma);
	
	do{
		Debug2("k: " << k << " sigma est: " << Sigma << " \npi est: " << pi.transpose() << " mu est: " << Mu.transpose());
		rice_estep(n, k, X, gamm, rgamm, pi, Mu, Sigma);
		Debug2("k: " << k << " sigma est: " << Sigma << " \npi est: " << pi.transpose() << " mu est: " << Mu.transpose());
		// show_head(gamm, 3);
		rice_mstep(X, n, k, pi, Mu, Sigma, gamm, rgamm);
		oldllhd = llhdval;
		llhdval = observedDataLogLikelihood(X, n, pi, Mu, k, Sigma);
		iter++;
		
	}
	while (((llhdval - oldllhd) > eps*std::fabs(oldllhd)) && (iter < maxiter));

	llhdval = observedDataLogLikelihood(X, n, pi, Mu, k, Sigma);
	BIC = 2*llhdval - (2*k)*std::log(n); 	// k + (k-1) + 1 // assuming last mu is 0
	class1 = classify(n, k, gamm);
	ICL = icl(n, k, X, class1, Mu, pi, Sigma);
	
	return class1;
}





// Change
double Est_var(Vector_eig r_col, int min_grp = 5, int max_grp = 15, int init_iter = 4, int EM_iter = 15, double eps = 0.0001){
	
	double sig_sq = 0.0;
	double sig = std::sqrt(sig_sq), best_sig;
	int J = 5, n = r_col.rows();
	double llhdval, ICL, BIC;
	Vector_eig BIC_vec(max_grp - min_grp + 1), Sig_vec(max_grp - min_grp + 1);
	
	// from 2 grp to max_grp
	for(J = 0; J < max_grp - min_grp + 1; ++J){
	
		Vector_eig pi = Vector_eig::Zero(J+min_grp), mu = Vector_eig::Zero(J+min_grp);
		//Est_var_init(r_col, J+min_grp, sig_sq, pi, mu, 5);
		//Mix_rice(r_col, J+min_grp, sig_sq, pi, mu, 5, 10);
		//random_inits(r_col, J+min_grp, pi, mu, sig);
		Est_var_init(r_col, J+min_grp, sig_sq, pi, mu, init_iter);
		sig = std::sqrt(sig_sq);
	
		rice_emcluster(n, J+min_grp, pi, r_col, mu, sig, EM_iter, eps, llhdval, ICL, BIC);
		Debug1( "\tJ+min_grp: " << J+min_grp << ", sigma est: " << sig << " \npi est: " << pi.transpose() << "\nmu est: " << mu.transpose() << "\n");
		BIC_vec(J) = BIC;
		Sig_vec(J) = sig;
	}
	
	
	
	//get location of minimum
	MatrixXd::Index maxRow, maxCol;
	double max = BIC_vec.maxCoeff(&maxRow, &maxCol);
	Debug1("BIC_vec: " << BIC_vec.transpose());
	Debug1("Sig_vec: " << Sig_vec.transpose());
	best_sig = Sig_vec(maxRow, maxCol);
	
	
	return best_sig;
}





int main(int argc, char * argv[]) {

	if (argc != 4) {
		fprintf(stderr, "\nUsage: %s <file_name> <SD_file_name> <will_write_to_a_file?> <temp_val> \n", argv[0]);
		exit(EXIT_FAILURE);
	}
	char *data_file, *sd_file;
	data_file = argv[1]; 	sd_file = argv[2]; 	char will_write = *(argv[3])-48;		// Converted from ascii
	short our_dim[8];
	
	


	
	Matrix_eig r = Preprocess_data(data_file, our_dim, will_write);
	Vector_eig sd_est = Vector_eig::Zero(r.cols());
	for(int i = 0; i < r.cols(); ++i){
		sd_est(i) = Est_var(r.col(i));
	}
	Debug0("sd_est: " << sd_est.transpose());
	
	std::ofstream myfile ("Dummy_sd.txt");
 	if (myfile.is_open()){
 		for(int i = 0; i < r.cols(); ++i){
 			myfile << sd_est(i) << "\n";
 		}
		myfile.close();
	}
	
	return 0;
}

























Eigen::VectorXi classify(int n, int k, Matrix_eig &Gamma){
	int i, j;
	double rowmax = -Inf;
	Eigen::VectorXi clas  = Eigen::VectorXi::Zero(n);
	for (i = 0; i < n; i++) {
		rowmax = -Inf;
		for (j = 0; j < k; j++) {
			if (Gamma(i, j) > rowmax) {
				clas[i] = j;
				rowmax = Gamma(i, j);
			}
		}
	}
	return clas;
}


double icl(int n, int k, const Vector_eig &X, Eigen::VectorXi &class1, 
			Vector_eig &Mu, Vector_eig &pi, double sigma){
	int i;
	double sum = 0.;
	for (i = 0; i < n; i++) {
		sum += std::log(pi[class1[i]]) - 2*std::log(sigma) 
			- (SQ(X[i]) - 2 * X[i] * Mu[class1[i]] + SQ(Mu[class1[i]]) )/(2*SQ(sigma));
	}
	return sum;
}


// from loglikelihood.c file
double observedDataLogLikelihood(const Vector_eig &y, const int numPoints,
								const Vector_eig &p, const Vector_eig &mu,
								const int numClusters,
								const double sigma){

  int i;
  double result = 0.0;
  const double sigmaSq = sigma * sigma;
  const double twoSigmaSqR = 1.0 / (2 * sigmaSq);

  for (i = 0; i < numPoints; ++i) {
    int j;
    const double yi = y[i];
    const double yiScaled = yi / sigmaSq;
    const double logYiScaled = std::log(yiScaled);
    const double ys = yi * yi * twoSigmaSqR;
    double sum = 0.0;

    /*
     * XXX! How can we make this safer numerically?
     */

    for (j = 0; j < numClusters; ++j) {
      const double muj = mu[j];
      sum += p[j] * std::exp(logYiScaled + logBesselI0(muj * yiScaled)
			- ys - muj * twoSigmaSqR * muj);
    }
    result += std::log(sum);
  }
  return result;
}

