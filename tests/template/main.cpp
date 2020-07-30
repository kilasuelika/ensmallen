/**
 * @file main.cpp
 * @author Zhou Yao
 *
 * ensmallen is free software; you may redistribute it and/or modify it under
 * the terms of the 3-clause BSD license.  You should have received a copy of
 * the 3-clause BSD license along with ensmallen.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */

#include<numeric>
#include<Eigen/Core>


#include<cmath>
#include<iostream>

#include "../../include/ensmallen_bits/eigen_helper/Rand.hpp"

#include "../../include/ensmallen_bits/gradient_descent/gradient_descent.hpp"
#include "../../include/ensmallen_bits/problems/ackley_function.hpp"


using namespace ens;
using namespace ens::test;
using namespace std;
using namespace Eigen;

int main(){
	AckleyFunction f;
	auto init=f.GetInitialPoint();
	
	GradientDescent opt;
	opt.Optimize(f,init);
	cout<<"Solution:"<<endl<<init<<endl;
	cout<<"Objective value: "<<f.Evaluate(init)<<endl;
	VectorXd g=init;
	f.Gradient(init,g);
	cout<<"Gradient: "<<g<<endl;
	
	VectorXd s(2);
	s<<0,0;
	cout<<"Target solution:"<<endl<<s<<endl;
	cout<<"Objective value: "<<f.Evaluate(s)<<endl;
}