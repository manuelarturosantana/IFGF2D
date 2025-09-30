#include "LinearOperator.hpp"

LinearOperator::LinearOperator(int rows, int cols, const std::function<void (Eigen::VectorXcd&, Eigen::VectorXcd&)>& mapping)
: rows_(rows), 
  cols_(cols),
  mapping_(mapping)
{
}

LinearOperator::~LinearOperator()
{
}

Eigen::VectorXcd LinearOperator::operator*(Eigen::VectorXcd &x) const
{

  Eigen::VectorXcd solution;

  mapping_(x, solution);

  return solution;

}
