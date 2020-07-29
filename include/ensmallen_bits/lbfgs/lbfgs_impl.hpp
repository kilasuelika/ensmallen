/**
 * @file lbfgs_impl.hpp
 * @author Dongryeol Lee (dongryel@cc.gatech.edu)
 * @author Ryan Curtin
 *
 * The implementation of the L_BFGS optimizer.
 *
 * ensmallen is free software; you may redistribute it and/or modify it under
 * the terms of the 3-clause BSD license.  You should have received a copy of
 * the 3-clause BSD license along with ensmallen.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef ENSMALLEN_LBFGS_LBFGS_IMPL_HPP
#define ENSMALLEN_LBFGS_LBFGS_IMPL_HPP

// In case it hasn't been included yet.
#include<vector>

namespace ens {

/**
 * Initialize the L_BFGS object.
 *
 * @param numBasis Number of memory points to be stored (default 5).
 * @param maxIterations Maximum number of iterations for the optimization
 *     (0 means no limit and may run indefinitely).
 * @param armijoConstant Controls the accuracy of the line search routine for
 *     determining the Armijo condition.
 * @param wolfe Parameter for detecting the Wolfe condition.
 * @param minGradientNorm Minimum gradient norm required to continue the
 *     optimization.
 * @param factr Minimum relative function value decrease to continue
 *     the optimization.
 * @param maxLineSearchTrials The maximum number of trials for the line search
 *     (before giving up).
 * @param minStep The minimum step of the line search.
 * @param maxStep The maximum step of the line search.
 */
inline L_BFGS::L_BFGS(const size_t numBasis,
                      const size_t maxIterations,
                      const double armijoConstant,
                      const double wolfe,
                      const double minGradientNorm,
                      const double factr,
                      const size_t maxLineSearchTrials,
                      const double minStep,
                      const double maxStep) :
    numBasis(numBasis),
    maxIterations(maxIterations),
    armijoConstant(armijoConstant),
    wolfe(wolfe),
    minGradientNorm(minGradientNorm),
    factr(factr),
    maxLineSearchTrials(maxLineSearchTrials),
    minStep(minStep),
    maxStep(maxStep),
    terminate(false)
{
  // Nothing to do.
}

/**
 * Calculate the scaling factor, gamma, which is used to scale the Hessian
 * approximation matrix.  See method M3 in Section 4 of Liu and Nocedal
 * (1989).
 *
 * @return The calculated scaling factor.
 * @param gradient The gradient at the initial point.
 * @param s Differences between the iterate and old iterate matrix.
 * @param y Differences between the gradient and the old gradient matrix.
 */
template<typename MatType>
double L_BFGS::ChooseScalingFactor(const size_t iterationNum,
                                   const MatType& gradient,
                                   const std::vector<MatType>& s,
                                   const std::vector<MatType>& y)
{
  typedef typename MatType::Scalar ElemType;

  double scalingFactor = 1.0;
  if (iterationNum > 0)
  {
    int previousPos = (iterationNum - 1) % numBasis;
    // Get s and y matrices once instead of multiple times.
	const auto& sMat=s[previousPos];
	const auto& yMat=y[previousPos];
	scalingFactor= sMat.dot(yMat)/yMat.dot(yMat);
  }
  else
  {
    //scalingFactor = 1.0 / sqrt(dot(gradient, gradient));
	scalingFactor = 1.0 / sqrt(gradient.dot(gradient));
  }

  return scalingFactor;
}

/**
 * Find the L_BFGS search direction.
 *
 * @param gradient The gradient at the current point.
 * @param iterationNum The iteration number.
 * @param scalingFactor Scaling factor to use (see ChooseScalingFactor_()).
 * @param s Differences between the iterate and old iterate matrix.
 * @param y Differences between the gradient and the old gradient matrix.
 * @param searchDirection Vector to store search direction in.
 */
template<typename MatType>
void L_BFGS::SearchDirection(const MatType& gradient,
                             const size_t iterationNum,
                             const double scalingFactor,
                             const std::vector<MatType>& s,
                             const std::vector<MatType>& y,
                             MatType& searchDirection)
{
  // Start from this point.
  searchDirection = gradient;

  // See "A Recursive Formula to Compute H * g" in "Updating quasi-Newton
  // matrices with limited storage" (Nocedal, 1980).
  typedef typename MatType::Scalar ElemType;

  // Temporary variables.
  //arma::Col<CubeElemType> rho(numBasis);
  //arma::Col<CubeElemType> alpha(numBasis);
	Eigen::Matrix<ElemType,Eigen::Dynamic,1> rho(numBasis), alpha(numBasis);
  size_t limit = (numBasis > iterationNum) ? 0 : (iterationNum - numBasis);
  for (size_t i = iterationNum; i != limit; i--)
  {
    int translatedPosition = (i + (numBasis - 1)) % numBasis;
    //rho.coeffRef(iterationNum - i) = 1.0 / arma::dot(y.slice(translatedPosition),
    //                                        s.slice(translatedPosition));
	rho.coeffRef(iterationNum - i) = 1.0 / y[translatedPosition].dot(s[translatedPosition]);
    //alpha[iterationNum - i] = rho[iterationNum - i] *
    //    arma::dot(s.slice(translatedPosition), searchDirection);
	alpha.coeffRef(iterationNum - i)=rho[iterationNum - i] *s[translatedPosition].dot(searchDirection);
    searchDirection -= alpha.coeff(iterationNum - i) * y[translatedPosition];
  }

  searchDirection *= scalingFactor;

  for (size_t i = limit; i < iterationNum; i++)
  {
    int translatedPosition = i % numBasis;
    double beta = rho.coeffRef(iterationNum - i - 1) *
        y[translatedPosition].dot(searchDirection);
    searchDirection += (alpha.coeff(iterationNum - i - 1) - beta) *
        s[translatedPosition];
  }

  // Negate the search direction so that it is a descent direction.
  searchDirection *= -1;
}

/**
 * Update the y and s matrices, which store the differences between
 * the iterate and old iterate and the differences between the gradient and the
 * old gradient, respectively.
 *
 * @param iterationNum Iteration number.
 * @param iterate Current point.
 * @param oldIterate Point at last iteration.
 * @param gradient Gradient at current point (iterate).
 * @param oldGradient Gradient at last iteration point (oldIterate).
 * @param s Differences between the iterate and old iterate matrix.
 * @param y Differences between the gradient and the old gradient matrix.
 */
template<typename MatType>
void L_BFGS::UpdateBasisSet(const size_t iterationNum,
                            const MatType& iterate,
                            const MatType& oldIterate,
                            const MatType& gradient,
                            const MatType& oldGradient,
                            std::vector<MatType>& s,
                            std::vector<MatType>& y)
{
  // Overwrite a certain position instead of pushing everything in the vector
  // back one position.
  int overwritePos = iterationNum % numBasis;
  s[overwritePos] = iterate - oldIterate;
  y[overwritePos] = gradient - oldGradient;
}

/**
 * Perform a back-tracking line search along the search direction to calculate a
 * step size satisfying the Wolfe conditions.
 *
 * @param function Function to optimize.
 * @param functionValue Value of the function at the initial point.
 * @param iterate The initial point to begin the line search from.
 * @param gradient The gradient at the initial point.
 * @param searchDirection A vector specifying the search direction.
 * @param finalStepSize The resulting step size used.
 * @param callbacks Callback functions.
 *
 * @return false if no step size is suitable, true otherwise.
 */
template<typename FunctionType,
         typename ElemType,
         typename MatType>
bool L_BFGS::LineSearch(FunctionType& function,
                        ElemType& functionValue,
                        MatType& iterate,
                        MatType& gradient,
                        MatType& newIterateTmp,
                        const MatType& searchDirection,
                        double& finalStepSize)
{
  // Default first step size of 1.0.
  double stepSize = 1.0;
  finalStepSize = 0.0; // Set only when we take the step.

  // The initial linear term approximation in the direction of the
  // search direction.
  ElemType initialSearchDirectionDotGradient =
      //arma::dot(gradient, searchDirection);
	  gradient.dot(searchDirection);

  // If it is not a descent direction, just report failure.
  if (initialSearchDirectionDotGradient > 0.0)
  {
    std::cout << "L-BFGS line search direction is not a descent direction "
        << "(terminating)!" << std::endl;
    return false;
  }

  // Save the initial function value.
  ElemType initialFunctionValue = functionValue;

  // Unit linear approximation to the decrease in function value.
  ElemType linearApproxFunctionValueDecrease = armijoConstant *
      initialSearchDirectionDotGradient;

  // The number of iteration in the search.
  size_t numIterations = 0;

  // Armijo step size scaling factor for increase and decrease.
  const double inc = 2.1;
  const double dec = 0.5;
  double width = 0;
  double bestStepSize = 1.0;
  ElemType bestObjective = std::numeric_limits<ElemType>::max();

  while (true)
  {
    // Perform a step and evaluate the gradient and the function values at that
    // point.
    newIterateTmp = iterate;
    newIterateTmp += stepSize * searchDirection;
    functionValue = function.EvaluateWithGradient(newIterateTmp, gradient);

    if (functionValue < bestObjective)
    {
      bestStepSize = stepSize;
      bestObjective = functionValue;
    }
    numIterations++;

    if (functionValue > initialFunctionValue + stepSize *
        linearApproxFunctionValueDecrease)
    {
      width = dec;
    }
    else
    {
      // Check Wolfe's condition.
      ElemType searchDirectionDotGradient = gradient.dot(searchDirection);

      if (searchDirectionDotGradient < wolfe *
          initialSearchDirectionDotGradient)
      {
        width = inc;
      }
      else
      {
        if (searchDirectionDotGradient > -wolfe *
            initialSearchDirectionDotGradient)
        {
          width = dec;
        }
        else
        {
          break;
        }
      }
    }

    // Terminate when the step size gets too small or too big or it
    // exceeds the max number of iterations.
    const bool cond1 = (stepSize < minStep);
    const bool cond2 = (stepSize > maxStep);
    const bool cond3 = (numIterations >= maxLineSearchTrials);
    if (cond1 || cond2 || cond3)
      break;

    // Scale the step size.
    stepSize *= width;
  }

  // Move to the new iterate.
  iterate += bestStepSize * searchDirection;
  finalStepSize = bestStepSize;
  return true;
}

/**
 * Use L_BFGS to optimize the given function, starting at the given iterate
 * point and performing no more than the specified number of maximum iterations.
 * The given starting point will be modified to store the finishing point of the
 * algorithm.
 *
 * @param numIterations Maximum number of iterations to perform
 * @param iterate Starting point (will be modified)
 * @param callbacks Callback functions.
 */
template<typename FunctionType,
         typename MatType>
typename MatType::Scalar L_BFGS::Optimize(FunctionType& function,
                 MatType& iterate)
{
	
	using namespace Eigen;
  // Convenience typedefs.
  typedef typename MatType::Scalar ElemType;
  
  // Ensure that the cubes holding past iterations' information are the right
  // size.  Also set the current best point value to the maximum.
  const size_t rows = iterate.rows();
  const size_t cols = iterate.cols();

  MatType newIterateTmp(rows, cols);
  //arma::Cube<ElemType> s(rows, cols, numBasis);
  std::vector<MatType> s(numBasis,MatType(rows)),y(numBasis,MatType(rows));

  // The old iterate to be saved.
  MatType oldIterate(iterate.rows());
  
  // Whether to optimize until convergence.
  bool optimizeUntilConvergence = (maxIterations == 0);

  // The gradient: the current and the old.
  MatType gradient(iterate.rows());
  MatType oldGradient(iterate.rows());

  // The search direction.
  MatType searchDirection(iterate.rows());

  // The initial function value and gradient.
  ElemType functionValue = function.EvaluateWithGradient(iterate, gradient);

  ElemType prevFunctionValue = functionValue;

  // The main optimization loop.
  for (size_t itNum = 0; (optimizeUntilConvergence || (itNum != maxIterations))
	  ; ++itNum)
  {
    prevFunctionValue = functionValue;

    // Break when the norm of the gradient becomes too small.
    //
    // But don't do this on the first iteration to ensure we always take at
    // least one descent step.
    if (itNum > 0 && gradient.norm() < minGradientNorm)
    {
      std::cout << "L-BFGS gradient norm too small (terminating successfully)."
          << std::endl;
      break;
    }

    // Break if the objective is not a number.
    if (std::isnan(functionValue))
    {
      std::cout << "L-BFGS terminated with objective " << functionValue << "; "
          << "are the objective and gradient functions implemented correctly?"
          << std::endl;
      break;
    }

    // Choose the scaling factor.
    double scalingFactor = ChooseScalingFactor(itNum, gradient, s, y);
    if (scalingFactor == 0.0)
    {
      std::cout << "L-BFGS scaling factor computed as 0 (terminating successfully)."
          << std::endl;
      break;
    }

    // Build an approximation to the Hessian and choose the search
    // direction for the current iteration.
    SearchDirection(gradient, itNum, scalingFactor, s, y, searchDirection);

    // Save the old iterate and the gradient before stepping.
    oldIterate = iterate;
    oldGradient = gradient;

    double stepSize; // Set by LineSearch().
    if (!LineSearch(function, functionValue, iterate, gradient, newIterateTmp,
        searchDirection, stepSize))
    {
      std::cout << "Line search failed.  Stopping optimization." << std::endl;
      break; // The line search failed; nothing else to try.
    }

    // It is possible that the difference between the two coordinates is zero.
    // In this case we terminate successfully.
    if (stepSize == 0.0)
    {
      std::cout << "L-BFGS step size of 0 (terminating successfully)."
          << std::endl;
      break;
    }

    // If we can't make progress on the gradient, then we'll also accept
    // a stable function value.
    const double denom = std::max(
        std::max(std::abs(prevFunctionValue), std::abs(functionValue)),
        (ElemType) 1.0);
    if ((prevFunctionValue - functionValue) / denom <= factr)
    {
      std::cout << "L-BFGS function value stable (terminating successfully)."
          << std::endl;
      break;
    }

    // Overwrite an old basis set.
    UpdateBasisSet(itNum, iterate, oldIterate, gradient, oldGradient, s, y);

    
  } // End of the optimization loop.

  return functionValue;
}

} // namespace ens

#endif // ENSMALLEN_LBFGS_LBFGS_IMPL_HPP

