/**
 * @file gradient_descent_impl.hpp
 * @author Sumedh Ghaisas
 * @author Zhou Yao
 * Simple gradient descent implementation.
 *
 * ensmallen is free software; you may redistribute it and/or modify it under
 * the terms of the 3-clause BSD license.  You should have received a copy of
 * the 3-clause BSD license along with ensmallen.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef ENSMALLEN_GRADIENT_DESCENT_GRADIENT_DESCENT_IMPL_HPP
#define ENSMALLEN_GRADIENT_DESCENT_GRADIENT_DESCENT_IMPL_HPP

// In case it hasn't been included yet.
#include "gradient_descent.hpp"

namespace ens {

//! Constructor.
inline GradientDescent::GradientDescent(
    const double stepSize,
    const size_t maxIterations,
    const double tolerance) :
    stepSize(stepSize),
    maxIterations(maxIterations),
    tolerance(tolerance)
{ /* Nothing to do. */ }

//! Optimize the function (minimize).
  template<typename FunctionType,
           typename MatType>
      typename MatType::Scalar
  GradientDescent::Optimize(FunctionType& function,
           MatType& iterate)
{
  // Convenience typedefs.
  typedef typename MatType::Scalar ElemType;

  // To keep track of where we are and how things are going.
  ElemType overallObjective = std::numeric_limits<ElemType>::max();
  ElemType lastObjective = std::numeric_limits<ElemType>::max();

  MatType gradient(iterate.rows());

  for (size_t i = 1; i != maxIterations; ++i)
  {
    overallObjective = function.EvaluateWithGradient(iterate, gradient);


    // Output current objective function.
    std::cout << "Gradient Descent: iteration " << i << ", objective "
        << overallObjective << "." << std::endl;

    if (std::isnan(overallObjective) || std::isinf(overallObjective))
    {
      std::cout << "Gradient Descent: converged to " << overallObjective
          << "; terminating" << " with failure.  Try a smaller step size?"
          << std::endl;

      return overallObjective;
    }

    if (std::abs(lastObjective - overallObjective) < tolerance)
    {
      std::cout << "Gradient Descent: minimized within tolerance "
          << tolerance << "; " << "terminating optimization." << std::endl;

      return overallObjective;
    }

    // Reset the counter variables.
    lastObjective = overallObjective;

    // And update the iterate.
    iterate -= stepSize * gradient;

  }

  std::cout << "Gradient Descent: maximum iterations (" << maxIterations
      << ") reached; " << "terminating optimization." << std::endl;

  return overallObjective;
}

  template<typename FunctionType,
           typename MatType>
 typename MatType::Scalar::type
  GradientDescent::Optimize(FunctionType& function,
           MatType& iterate,
           const std::vector<bool>& categoricalDimensions,
           const Eigen::VectorXi& numCategories)
{
  if (categoricalDimensions.size() != iterate.rows())
  {
    std::ostringstream oss;
    oss << "GradientDescent::Optimize(): expected information about "
        << iterate.n_rows << " dimensions in categoricalDimensions, "
        << "but got " << categoricalDimensions.size();
    throw std::invalid_argument(oss.str());
  }

  if (numCategories.rows() != iterate.rows())
  {
    std::ostringstream oss;
    oss << "GradientDescent::Optimize(): expected numCategories to have length "
        << "equal to number of dimensions (" << iterate.n_rows << ") but it has"
        << " length " << numCategories.rows();
    throw std::invalid_argument(oss.str());
  }

  for (size_t i = 0; i < categoricalDimensions.size(); ++i)
  {
    if (categoricalDimensions[i])
    {
      std::ostringstream oss;
      oss << "GradientDescent::Optimize(): the dimension " << i
          << "is not numeric" << std::endl;
      throw std::invalid_argument(oss.str());
    }
  }

  return Optimize(function, iterate);
}

} // namespace ens

#endif
