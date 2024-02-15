#include "heat_geodesics.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>
#include <igl/isolines_map.h>
#include <igl/edge_lengths.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
  // Load input meshes
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  double lambda = 1e-5;
  igl::read_triangle_mesh(
    (argc>1?argv[1]:"../data/penguin.obj"),V,F);

  // Some function to upsample
  igl::opengl::glfw::Viewer viewer;
  std::cout<<R"(
  S,s  select a new source point
  H,h  heat geodesic distances
  L    lighting
)";
  viewer.data().show_lines = false;
  viewer.data().set_mesh(V,F);
  Eigen::VectorXd P(V.rows());
  P.setZero();
  const auto & update = [&]()
  {
    Eigen::MatrixXd CM;
    igl::parula(Eigen::VectorXd::LinSpaced(50,0,1).eval(),false,CM);
    igl::isolines_map(Eigen::MatrixXd(CM),CM);
    viewer.data().set_colormap(CM);
    viewer.data().set_data(P);
  };
  int s = 0;
  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int)
  {
    Eigen::VectorXi J, I;
    switch(key)
    {
      case 'S':
      case 's':
        s = rand() % V.rows();
        std::cout << "New source point: " << s << std::endl;
        viewer.data().set_points(V.row(s),
          Eigen::RowVector3d(8./255.,81./255.,156./255.));
        break;
      case 'H':
      case 'h':
        //////////////////////////////////////////////////////////////////////
        // Heat geodesics
        //////////////////////////////////////////////////////////////////////
        {
          Eigen::MatrixXd lE;
          igl::edge_lengths(V, F, lE);
          const double ll = lE.mean();
          Eigen::VectorXi gamma(1);
          gamma(0) = s;
          heat_geodesics(V, F, gamma, ll*ll, P);
        }
        break;
      case 'L':
        // Toggle lighting
        viewer.core().lighting_factor = 1.0- viewer.core().lighting_factor;
        break;
      default:
        return false;
    }
    update();
    return true;
  };

  update();
  viewer.callback_key_pressed(viewer, 's', 0);
  viewer.launch();

  return EXIT_SUCCESS;
}
