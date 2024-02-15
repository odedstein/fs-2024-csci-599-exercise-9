#ifndef HEAT_GEODESICS_H
#define HEAT_GEODESICS_H
#include <Eigen/Core>
  /// Compute fast approximate geodesic distances from a
  /// set of selected source vertices (gamma).
  ///
  /// @param[in] V  #V by dim list of mesh vertex positions
  /// @param[in] F  #F by 3 list of mesh face indices into V
  /// @param[in] gamma  #gamma list of indices into V of source vertices
  /// @param[in] t  "heat" parameter (smaller --> more accurate, less stable)
  /// @param[out] D  #V list of distances to gamma 
  ///
  /// \fileinfo
  void heat_geodesics(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::VectorXi & gamma,
    const double t,
    Eigen::VectorXd & D);

#endif
