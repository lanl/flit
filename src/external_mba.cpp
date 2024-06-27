//
// © 2024. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 
// for Los Alamos National Laboratory (LANL), which is operated by 
// Triad National Security, LLC for the U.S. Department of Energy/National Nuclear 
// Security Administration. All rights in the program are reserved by 
// Triad National Security, LLC, and the U.S. Department of Energy/National 
// Nuclear Security Administration. The Government is granted for itself and 
// others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide 
// license in this material to reproduce, prepare. derivative works, 
// distribute copies to the public, perform publicly and display publicly, 
// and to permit others to do so.
//
// Author:
//    Kai Gao, kaigao@lanl.gov
//

#include <mba.hpp>

using namespace std;

extern "C" {

    void mba1(int n,
              float *x,
              float *y,
              int nn,
              float *xx,
              float *yy);

    void mba2(int n,
              float *x,
              float *y,
              float *z,
              int nn,
              float *xx,
              float *yy,
              float *zz);

    void mba3(int n,
              float *x,
              float *y,
              float *z,
              float *v,
              int nn,
              float *xx,
              float *yy,
              float *zz,
              float *vv);

    void mba1_double(int n,
                     double *x,
                     double *y,
                     int nn,
                     double *xx,
                     double *yy);

    void mba2_double(int n,
                     double *x,
                     double *y,
                     double *z,
                     int nn,
                     double *xx,
                     double *yy,
                     double *zz);

    void mba3_double(int n,
                     double *x,
                     double *y,
                     double *z,
                     double *v,
                     int nn,
                     double *xx,
                     double *yy,
                     double *zz,
                     double *vv);

}

void mba1(int n,
          float *x,
          float *y,
          int nn,
          float *xx,
          float *yy)
{

    // Coordinates of the scattered data points
    vector<mba::point<1>> coo;
    vector<double> val;

    for (size_t i = 0; i < n; i++)
    {
        coo.push_back(mba::point<1> {(double)x[i]});
        val.push_back((double)y[i]);
    }

    // Bounding box containing the data points, must > real range
    double xmin = std::min(*std::min_element(x, x + n), *std::min_element(xx, xx + nn));
    double xmax = std::max(*std::max_element(x, x + n), *std::max_element(xx, xx + nn));

    mba::point<1> lo = {xmin - 0.05 * (xmax - xmin)};
    mba::point<1> hi = {xmax + 0.05 * (xmax - xmin)};

    // Initial grid size.
    mba::index<1> grid = {3};

    // Algorithm setup.
    mba::MBA<1> interp(lo, hi, grid, coo, val);

    // Get interpolated value at arbitrary location.
    #pragma omp parallel for
    for (size_t i = 0; i < nn; i++)
    {
        yy[i] = float(interp(mba::point<1> {xx[i]}));
    }

    coo.clear();
    val.clear();
}

void mba2(int n,
          float *x,
          float *y,
          float *z,
          int nn,
          float *xx,
          float *yy,
          float *zz)
{

    // Coordinates of the scattered data points
    vector<mba::point<2>> coo;
    vector<double> val;

    for (size_t i = 0; i < n; i++)
    {
        coo.push_back(mba::point<2> {(double)x[i], (double)y[i]});
        val.push_back((double)z[i]);
    }

    // Bounding box containing the data points, must > real range
    double xmin = std::min(*std::min_element(x, x + n), *std::min_element(xx, xx + nn));
    double xmax = std::max(*std::max_element(x, x + n), *std::max_element(xx, xx + nn));
    double ymin = std::min(*std::min_element(y, y + n), *std::min_element(yy, yy + nn));
    double ymax = std::max(*std::max_element(y, y + n), *std::max_element(yy, yy + nn));

    mba::point<2> lo = {xmin - 0.05 * (xmax - xmin),
                        ymin - 0.05 * (ymax - ymin)
                       };
    mba::point<2> hi = {xmax + 0.05 * (xmax - xmin),
                        ymax + 0.05 * (ymax - ymin)
                       };

    // Initial grid size.
    mba::index<2> grid = {3, 3};

    // Algorithm setup.
    mba::MBA<2> interp(lo, hi, grid, coo, val);

    // Get interpolated value at arbitrary location.
    #pragma omp parallel for
    for (size_t i = 0; i < nn; i++)
    {
        zz[i] = float(interp(mba::point<2> {xx[i], yy[i]}));
    }

    coo.clear();
    val.clear();
}

void mba3(int n,
          float *x,
          float *y,
          float *z,
          float *v,
          int nn,
          float *xx,
          float *yy,
          float *zz,
          float *vv)
{

    // Coordinates of the scattered data points
    vector<mba::point<3>> coo;
    vector<double> val;

    for (size_t i = 0; i < n; i++)
    {
        coo.push_back(mba::point<3> {(double)x[i], (double)y[i], (double)z[i]});
        val.push_back((double)z[i]);
    }

    // Bounding box containing the data points, must > real range
    double xmin = std::min(*std::min_element(x, x + n), *std::min_element(xx, xx + nn));
    double xmax = std::max(*std::max_element(x, x + n), *std::max_element(xx, xx + nn));
    double ymin = std::min(*std::min_element(y, y + n), *std::min_element(yy, yy + nn));
    double ymax = std::max(*std::max_element(y, y + n), *std::max_element(yy, yy + nn));
    double zmin = std::min(*std::min_element(z, z + n), *std::min_element(zz, zz + nn));
    double zmax = std::max(*std::max_element(z, z + n), *std::max_element(zz, zz + nn));

    mba::point<3> lo = {xmin - 0.05 * (xmax - xmin),
                        ymin - 0.05 * (ymax - ymin),
                        zmin - 0.05 * (zmax - zmin)
                       };
    mba::point<3> hi = {xmax + 0.05 * (xmax - xmin),
                        ymax + 0.05 * (ymax - ymin),
                        zmax + 0.05 * (zmax - zmin)
                       };

    // Initial grid size.
    mba::index<3> grid = {3, 3, 3};

    // Algorithm setup.
    mba::MBA<3> interp(lo, hi, grid, coo, val);

    // Get interpolated value at arbitrary location.
    #pragma omp parallel for
    for (size_t i = 0; i < nn; i++)
    {
        vv[i] = float(interp(mba::point<3> {xx[i], yy[i], zz[i]}));
    }

    coo.clear();
    val.clear();
}

void mba1_double(int n,
                 double *x,
                 double *y,
                 int nn,
                 double *xx,
                 double *yy)
{

    // Coordinates of the scattered data points
    vector<mba::point<1>> coo;
    vector<double> val;

    for (size_t i = 0; i < n; i++)
    {
        coo.push_back(mba::point<1> {x[i]});
        val.push_back(y[i]);
    }

    // Bounding box containing the data points, must > real range
    double xmin = std::min(*std::min_element(x, x + n), *std::min_element(xx, xx + nn));
    double xmax = std::max(*std::max_element(x, x + n), *std::max_element(xx, xx + nn));

    mba::point<1> lo = {xmin - 0.05 * (xmax - xmin)};
    mba::point<1> hi = {xmax + 0.05 * (xmax - xmin)};

    // Initial grid size.
    mba::index<1> grid = {3};

    // Algorithm setup.
    mba::MBA<1> interp(lo, hi, grid, coo, val);

    // Get interpolated value at arbitrary location.
    #pragma omp parallel for
    for (size_t i = 0; i < nn; i++)
    {
        yy[i] = interp(mba::point<1> {xx[i]});
    }

    coo.clear();
    val.clear();
}

void mba2_double(int n,
                 double *x,
                 double *y,
                 double *z,
                 int nn,
                 double *xx,
                 double *yy,
                 double *zz)
{

    // Coordinates of the scattered data points
    vector<mba::point<2>> coo;
    vector<double> val;

    for (size_t i = 0; i < n; i++)
    {
        coo.push_back(mba::point<2> {x[i], y[i]});
        val.push_back(z[i]);
    }

    // Bounding box containing the data points, must > real range
    double xmin = std::min(*std::min_element(x, x + n), *std::min_element(xx, xx + nn));
    double xmax = std::max(*std::max_element(x, x + n), *std::max_element(xx, xx + nn));
    double ymin = std::min(*std::min_element(y, y + n), *std::min_element(yy, yy + nn));
    double ymax = std::max(*std::max_element(y, y + n), *std::max_element(yy, yy + nn));

    mba::point<2> lo = {xmin - 0.05 * (xmax - xmin),
                        ymin - 0.05 * (ymax - ymin)
                       };
    mba::point<2> hi = {xmax + 0.05 * (xmax - xmin),
                        ymax + 0.05 * (ymax - ymin)
                       };

    // Initial grid size.
    mba::index<2> grid = {3, 3};

    // Algorithm setup.
    mba::MBA<2> interp(lo, hi, grid, coo, val);

    // Get interpolated value at arbitrary location.
    #pragma omp parallel for
    for (size_t i = 0; i < nn; i++)
    {
        zz[i] = interp(mba::point<2> {xx[i], yy[i]});
    }

    coo.clear();
    val.clear();
}

void mba3_double(int n,
                 double *x,
                 double *y,
                 double *z,
                 double *v,
                 int nn,
                 double *xx,
                 double *yy,
                 double *zz,
                 double *vv)
{

    // Coordinates of the scattered data points
    vector<mba::point<3>> coo;
    vector<double> val;

    for (size_t i = 0; i < n; i++)
    {
        coo.push_back(mba::point<3> {x[i], y[i], z[i]});
        val.push_back(z[i]);
    }

    // Bounding box containing the data points, must > real range
    double xmin = std::min(*std::min_element(x, x + n), *std::min_element(xx, xx + nn));
    double xmax = std::max(*std::max_element(x, x + n), *std::max_element(xx, xx + nn));
    double ymin = std::min(*std::min_element(y, y + n), *std::min_element(yy, yy + nn));
    double ymax = std::max(*std::max_element(y, y + n), *std::max_element(yy, yy + nn));
    double zmin = std::min(*std::min_element(z, z + n), *std::min_element(zz, zz + nn));
    double zmax = std::max(*std::max_element(z, z + n), *std::max_element(zz, zz + nn));

    mba::point<3> lo = {xmin - 0.05 * (xmax - xmin),
                        ymin - 0.05 * (ymax - ymin),
                        zmin - 0.05 * (zmax - zmin)
                       };
    mba::point<3> hi = {xmax + 0.05 * (xmax - xmin),
                        ymax + 0.05 * (ymax - ymin),
                        zmax + 0.05 * (zmax - zmin)
                       };

    // Initial grid size.
    mba::index<3> grid = {3, 3, 3};

    // Algorithm setup.
    mba::MBA<3> interp(lo, hi, grid, coo, val);

    // Get interpolated value at arbitrary location.
    #pragma omp parallel for
    for (size_t i = 0; i < nn; i++)
    {
        vv[i] = interp(mba::point<3> {xx[i], yy[i], zz[i]});
    }

    coo.clear();
    val.clear();
}
