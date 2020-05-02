/**
 *           c++11-only implementation of the L-BFGS-B algorithm
 *
 * Copyright (c) 2014 Patrick Wieschollek
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <stdexcept>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <vector>

struct Options{
	double tol;
	double functol;
	double constrtol;
	int maxIter;
	int m ;
	Options(){
		tol = 1e-8;			//Subrata 1e-4
		functol = 1e-8;
		constrtol = 1e-2;
		maxIter = 1e4;
		m = 10;
	}
};

const double EPS=2.2204e-016;

typedef Eigen::MatrixXf Matrix_eig;
typedef Eigen::VectorXf Vector_eig;
typedef Eigen::VectorXf::Scalar Scalar_eig;

typedef Eigen::Map<Matrix_eig> ViewMatrix;
typedef Eigen::Map<Vector_eig> ViewVector;

typedef const Eigen::Map<const Matrix_eig> ReadMatrix;
typedef const Eigen::Map<const Vector_eig> ReadVector;

typedef std::function<double(const Eigen::VectorXf &x)> FunctionOracleType;
typedef std::function<void(const Eigen::VectorXf &x, Eigen::VectorXf &gradient)> GradientOracleType;

#define INF HUGE_VAL
#define Assert(x,m) if (!(x)) { throw (std::runtime_error(m)); }

//#define FAST

#ifdef FAST
#define Debug(x)
#else
//#define Debug(x) if(false){std::cout << "DEBUG: "<< x;std::cout<< std::endl;}
#define Debug(x) {std::cout << "DEBUG: "<< x; std::cout<< std::endl;}		//Subrata
#endif


#define min_eig(a,b) (((a)<(b))?(a):(b))
#define max_eig(a,b) (((a)>(b))?(a):(b))

std::vector<int> sort_indexes(const std::vector< std::pair<int,double> > &v) {
  std::vector<int> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i)
	  idx[i] = v[i].first;
  sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1].second < v[i2].second;});
  return idx;
}


#endif /* DEFINITIONS_H_ */
