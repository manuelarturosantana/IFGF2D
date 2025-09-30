#pragma once
// Make compatable with Eigen GMRES

#include <complex>
#include <functional>
#include <vector>
#include <Eigen/Geometry>

// Forward declare the operator class
class LinearOperator;

//--- Begin Eigen GMRES Compatability ------------------------------------------
namespace Eigen {
  namespace internal {

    // Eigen expects this type traits class to exist for use internally
    template<class _MatrixType>
    struct traits
    {
      typedef typename Eigen::Dense StorageKind;
    };

  }
}
//--- End Eigen GMRES Compatability --------------------------------------------

class LinearOperator : public Eigen::EigenBase<LinearOperator>
{
//--- Begin Eigen GMRES Compatability ------------------------------------------
public:
  typedef std::complex<double> Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  typedef Eigen::Dense StorageKind;

  enum { ColsAtCompileTime = 0 };
  enum { MaxColsAtCompileTime = 0 };

  template<typename Dest>
  void evalTo(Dest& dst) const { };

  // Needed for GMRES to determine operator dimensions
  int rows() const { return rows_; };
  int cols() const { return cols_; };

protected:
  int rows_;
  int cols_;
  std::function<void (const Eigen::VectorXcd&, Eigen::VectorXcd&)> mapping_;

//--- End Eigen GMRES Compatability --------------------------------------------

public:
  LinearOperator(int rows, int cols, const std::function<void (const Eigen::VectorXcd&, Eigen::VectorXcd&)>& mapping);
  
  ~LinearOperator();

  Eigen::VectorXcd operator*(const Eigen::VectorXcd &x) const;

};