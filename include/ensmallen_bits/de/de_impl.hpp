/**
 * @file de_impl.hpp
 * @author Rahul Ganesh Prabhu
 * @author Zhou Yao
 *
 * Implementation of Differential Evolution an evolutionary algorithm used for
 * global optimization of arbitrary functions.
 *
 * ensmallen is free software; you may redistribute it and/or modify it under
 * the terms of the 3-clause BSD license.  You should have received a copy of
 * the 3-clause BSD license along with ensmallen.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef ENSMALLEN_DE_DE_IMPL_HPP
#define ENSMALLEN_DE_DE_IMPL_HPP

#include "de.hpp"
#include<random>

namespace ens {

inline DE::DE(const size_t populationSize ,
              const size_t maxGenerations,
              const double crossoverRate,
              const double differentialWeight,
              const double tolerance):
    populationSize(populationSize),
    maxGenerations(maxGenerations),
    crossoverRate(crossoverRate),
    differentialWeight(differentialWeight),
    tolerance(tolerance)
{ /* Nothing to do here. */ }

//!Optimize the function
template<typename FunctionType,
         typename MatType>
typename MatType::Scalar DE::Optimize(FunctionType& function,
                                         MatType& iterateIn)
{
	using namespace Eigen;
	const int nargs=iterateIn.rows()*iterateIn.cols();

  // Convenience typedefs.
  typedef typename MatType::Scalar ElemType;
  
  // Population matrix. Each column is a candidate.
  std::vector<MatType> population;
  population.resize(populationSize);
  
  // Vector of fitness values corresponding to each candidate.
  Matrix<ElemType, Dynamic,1> fitnessValues(populationSize);
  
  // Population Size must be at least 3 for DE to work.
  if (populationSize < 3)
  {
    throw std::logic_error("CNE::Optimize(): population size should be at least"
        " 3!");
  }

  // Initialize helper variables.
  ElemType lastBestFitness = std::numeric_limits<ElemType>::max();
  MatType bestElement;

  // Generate a population based on a Gaussian distribution around the given
  // starting point. Also finds the best element of the population.
  eigen_helper::RandGen gen;
  for (size_t i = 0; i < populationSize; i++)
  {
	population[i].resize(nargs);
	gen.randn(population[i]);
    population[i] += iterateIn;
    fitnessValues.coeffRef(i) = function.Evaluate(population[i]);

    if (fitnessValues.coeff(i) < lastBestFitness)
    {
      lastBestFitness = fitnessValues.coeff(i);
      bestElement = population[i];
    }
  }

  // Iterate until maximum number of generations are completed.
  //For random integer.
  std::random_device rd;
  std::uniform_int_distribution<> dist(0, populationSize-1);
    
  for (size_t gen = 0; gen < maxGenerations; gen++)
  {
	  std::cout<<gen<<std::endl;
    // Generate new population based on /best/1/bin strategy.
    for (size_t member = 0; member < populationSize; member++)
    {
      iterateIn = population[member];

      // Generate two different random numbers to choose two random members.
      size_t l = 0, m = 0;
      do
      {
        l=dist(rd);
      }
      while (l == member);

      do
      {
		m=dist(rd);
      }
      while (m == member && m == l);

      // Generate new "mutant" from two randomly chosen members.
      MatType mutant = bestElement + differentialWeight *
          (population[l] - population[m]);

      // Perform crossover.
	  VectorXd cr;
	  cr.setRandom(nargs);
      for (size_t it = 0; it < iterateIn.rows(); it++)
      {
        if (cr.coeff(it) >= crossoverRate)
        {
          mutant.coeffRef(it) = iterateIn.coeff(it);
        }
      }

      ElemType iterateValue = function.Evaluate(iterateIn);
      
      const ElemType mutantValue = function.Evaluate(mutant);
      
      // Replace the current member if mutant is better.
      if (mutantValue < iterateValue)
      {
        iterateIn = mutant;
        iterateValue = mutantValue;
   
      }

      fitnessValues.coeffRef(member) = iterateValue;
      population[member] = iterateIn;
    }

    // Check for termination criteria.
    if (std::abs(lastBestFitness - fitnessValues.minCoeff()) < tolerance)
    {
      std::cout << "DE: minimized within tolerance " << tolerance << "; "
          << "terminating optimization." << std::endl;
      break;
    }

    // Update helper variables.
	int k;
	lastBestFitness = fitnessValues.minCoeff(&k);
	bestElement=population[k];
	
  }

  iterateIn = bestElement;

  return lastBestFitness;
}

} // namespace ens

#endif
