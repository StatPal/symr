/*
 * Compile and Run:
 g++ simple_subrata_bdd_2.cpp -o test -I /usr/include/eigen3 -O2 && ./test
 *
 *
 */


#include <iostream>
#include <iomanip>
#include "include/cppoptlib/meta.h"
#include "include/cppoptlib/boundedproblem.h"
#include "include/cppoptlib/solver/lbfgsbsolver.h"


// nolintnextline
//using namespace cppoptlib;		--Subrata commented out


/*Subrata - defined this function for testing const works or not in x - Works I guess*/
double temp_fn(const Eigen::VectorXd &x, const Eigen::VectorXd &r, double r1){
	return (5*x[0]*x[0] + 10*x[1]*x[1] + 50*sin(x[2]-1.570796) + r[0]*x[1] + r1 * x[0]);
}


template<typename T>
class Simple_Subrata_BDD : public cppoptlib::BoundedProblem<T> {
  public:
	using typename cppoptlib::BoundedProblem<T>::TVector;
	//const TVector r;
	TVector r;

  public:
	Simple_Subrata_BDD(const TVector y_) :											// Is it constructor?
		cppoptlib::BoundedProblem<T>(y_.size()), r(y_){}							// ' r(y_){}' was in next line.

	T r1;
	
	T value(const TVector &x) {												// Basically Eigen::VectorXd/VectorXf according to T
		double temp = 5*x[0]*x[0] + 10*x[1]*x[1] + 50*sin(x[2]-1.570796) + r[0]*x[1] + r1 * x[0];
		//double temp = temp_fn(x, r, r1);  // Subrata - check
		
		std::cout << "value:" << temp << "\n";
		return temp;
		// 5*x_0^2 + 100*x_1^2+ 50*sin(x_2-pi/2) + r*x_1 + r1*x_0
	}

	void gradient(const TVector &x, TVector &grad) {
		grad[0]  = 2*5*x[0] + r1;
		grad[1]  = 2*10*x[1] + r[0];
		grad[2]  = 50*cos(x[2]-1.570796);
	}
	
	// Possibly, we should add a destructor.
};




int op(Eigen::VectorXd x){
	for(int j = 0; j < 1; ++j){
		Eigen::VectorXd r = Eigen::VectorXd::Ones(4) * j;
	
		Simple_Subrata_BDD<double> f(r);
		Eigen::VectorXd lb(3), ub(3); 
		lb << 0.25, 0, 0;
		ub << 5, 1, 1;
		f.setLowerBound(lb);
		f.setUpperBound(ub);
		f.r1 = 0;
		
		//Print initial values:
		std::cout << "x at first: " << x.transpose() << ";\tf(x) at first: " << f(x) <<  std::endl;
		
		//Solve:
		cppoptlib::LbfgsbSolver<Simple_Subrata_BDD<double>> solver;
		
		
		cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
		// Create a Criteria class to set the solver's stop conditions
		crit.iterations = 500;														// Change the number of allowed iterations
		solver.setStopCriteria(crit);
		
		
		std::cout << "(" << std::setw(2) << crit.iterations << ")" 
			<< " ||dx|| = " << std::fixed << std::setw(8) << std::setprecision(4) << crit.gradNorm 
			<< std::endl;
		
		
		
		// Check Derivative:
		bool probably_correct_2 = f.checkGradient(x);
		if(probably_correct_2){
			std::cout << "\nTrue in Deriv" << "\n\n";
		} else {
			std::cout << "\nFalse in Deriv" << "\n\n";
		}
		
		
		solver.minimize(f, x);
		std::cout << "argmin: " << x.transpose() << ";\tf(x) in argmin: " << f(x) << std::endl;
		Eigen::VectorXd temp_grad = x; 		f.gradient(x, temp_grad);
		std::cout << "Final grad:" << temp_grad.transpose() << "\n";
		std::cout << "Solver status: " << solver.status() << std::endl;	//Guess: bad reports: under constraints => grad is not ~0 
		std::cout << "Final criteria values: " << std::endl << solver.criteria() << std::endl;
		
		
		std::cout << "Change values of parameters." << "\n";
		r = Eigen::VectorXd::Ones(4) * (j-1);
		f.r = r;
		f.r1 = -1;
		solver.minimize(f, x);
		std::cout << "After changing r and r1: " << "\n";
		std::cout << "argmin: " << x.transpose() << std::endl;
		std::cout << "f(x) in argmin: " << f(x) << std::endl;
		std::cout << "Solver status: " << solver.status() << std::endl;	//Guess: bad reports: under constraints => grad is not ~0 
		std::cout << "Final criteria values: " << std::endl << solver.criteria() << std::endl;
		
	}
	return 0;
}


int main() {
	
	Eigen::VectorXd x(3);
	x <<  1, 2, 2;
	x << 1, 0.5, 0.5;
	std::cin >> x[0];
	
	op(x);
	return 0;
}
