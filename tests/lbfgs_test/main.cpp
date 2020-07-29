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

#include "../../include/ensmallen_bits/lbfgs/lbfgs.hpp"
#include "../../include/ensmallen_bits/problems/ackley_function.hpp"


using namespace ens;
using namespace ens::test;
using namespace std;

int main(){
	AckleyFunction f;
	auto init=f.GetInitialPoint();
	init(0)=-1;init(1)=1;
	L_BFGS opt;
	opt.Optimize(f,init);
	//Eigen::VectorXd g(2);
	//f.EvaluateWithGradient(init,g);
	cout<<"Solution:"<<endl<<init<<endl;
	cout<<"Objective value: "<<f.Evaluate(init)<<endl;
	//f.EvaluateWithGradient(init,g);
	//cout<<"Gradient: "<<endl<<g<<endl;

	Eigen::VectorXd s(2);
	s<<0,0;
	cout<<"Target solution:"<<endl<<s<<endl;
	cout<<"Objective value: "<<f.Evaluate(s)<<endl;
}
