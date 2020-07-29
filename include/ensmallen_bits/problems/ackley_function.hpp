/**
 * @file ackley_function.hpp
 * @author Suryoday Basak
 *
 * Definition of the Ackley function.
 *
 * ensmallen is free software; you may redistribute it and/or modify it under
 * the terms of the 3-clause BSD license.  You should have received a copy of
 * the 3-clause BSD license along with ensmallen.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef ENSMALLEN_PROBLEMS_ACKLEY_FUNCTION_HPP
#define ENSMALLEN_PROBLEMS_ACKLEY_FUNCTION_HPP

#define _USE_MATH_DEFINES
#include <cmath>

namespace ens {
namespace test {

/**
 * The Ackley function, defined by
 *
 * \f[
 * f(x_1,x_2) = -20 * e^(-0.2 * sqrt(0.5 * (x_1^2 + x_2^2))) -
 *    e * (0.5(cos(2 * pi * x_1) + cos(2 * pi * x_2))) + e + 20
 * \f]
 *
 * This should optimize to f(x) = 0, at x = [0, 0].
 *
 * For more information, please refer to:
 *
 * @code
 * @book{Ackley1987,
 *   doi       = {10.1007/978-1-4613-1997-9},
 *   url       = {https://doi.org/10.1007/978-1-4613-1997-9},
 *   year      = {1987},
 *   publisher = {Springer {US}},
 *   author    = {David H. Ackley},
 *   title     = {A Connectionist Machine for Genetic Hillclimbing}
 * }
 * @endcode
 */
class AckleyFunction
{
 public:
  /**
   * Initialize the AckleyFunction.
   *
   * @param c Multiplicative constant with a default value of 2 * pi.
   * @param epsilon Coefficient to avoid division by zero (numerical stability).
   */
  AckleyFunction(const double c = 2 * M_PI,
                 const double epsilon = 1e-8);

  /**
   * Shuffle the order of function visitation. This may be called by the
   * optimizer.
   */
  void Shuffle();

  //! Return 1 (the number of functions).
  size_t NumFunctions() const { return 1; }

  //! Get the starting point.
  template<typename MatType = Eigen::VectorXd>
  MatType GetInitialPoint() const {
	  MatType m(2);
	  m<<-5,5;
	return m; 
	}

  /**
   * Evaluate a function for a particular batch-size.
   *
   * @param coordinates The function coordinates.
   * @param begin The first function.
   * @param batchSize Number of points to process.
   */
  template<typename MatType>
  typename MatType::Scalar Evaluate(const MatType& coordinates,
                                       const size_t begin,
                                       const size_t batchSize) const;

  /**
   * Evaluate a function with the given coordinates.
   *
   * @param coordinates The function coordinates.
   */
  template<typename MatType>
  typename MatType::Scalar Evaluate(const MatType& coordinates) const;

  /**
   * Evaluate the gradient of a function for a particular batch-size.
   *
   * @param coordinates The function coordinates.
   * @param begin The first function.
   * @param gradient The function gradient.
   * @param batchSize Number of points to process.
   */
  template<typename MatType>
  void Gradient(const MatType& coordinates,
                const size_t begin,
                MatType& gradient,
                const size_t batchSize) const;

  /**
   * Evaluate the gradient of a function with the given coordinates.
   *
   * @param coordinates The function coordinates.
   * @param gradient The function gradient.
   */
  template<typename MatType>
  void Gradient(const MatType& coordinates, MatType& gradient);

template<typename MatType>
 typename MatType::Scalar EvaluateWithGradient(const MatType& coordinates, MatType& gradient);

  //! Get the value used for c.
  double MultiplicativeConstant() const { return c; }
  //! Modify the value used for c.
  double& MultiplicativeConstant() { return c; }

  //! Get the value used for numerical stability.
  double Epsilon() const { return epsilon; }
  //! Modify the value used for numerical stability.
  double& Epsilon() { return epsilon; }

 private:
  //! The value of the multiplicative constant.
  double c;
  //! The value used for numerical stability.
  double epsilon;
};

} // namespace test
} // namespace ens

// Include implementation.
#include "ackley_function_impl.hpp"

#endif // ENSMALLEN_PROBLEMS_ACKLEY_FUNCTION_HPP
