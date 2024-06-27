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

#include "Splines.hh"

// external interface
extern "C" {

    void lspline1d(int n,
                   float *x,
                   float *y,
                   int nn,
                   float *xx,
                   float *yy);

    void cspline1d(int n,
                   float *x,
                   float *y,
                   int nn,
                   float *xx,
                   float *yy);

    void deriv_cspline1d(int n,
                         float *x,
                         float *y,
                         int nn,
                         float *xx,
                         float *yy,
                         int deriv);

    void deriv_pchip1d(int n,
                       float *x,
                       float *y,
                       int nn,
                       float *xx,
                       float *yy,
                       int deriv);

    void pchip1d(int n,
                 float *x,
                 float *y,
                 int nn,
                 float *xx,
                 float *yy);

    void quintic1d(int n,
                   float *x,
                   float *y,
                   int nn,
                   float *xx,
                   float *yy);

    void lspline1d_double(int n,
                          double *x,
                          double *y,
                          int nn,
                          double *xx,
                          double *yy);

    void cspline1d_double(int n,
                          double *x,
                          double *y,
                          int nn,
                          double *xx,
                          double *yy);

    void deriv_cspline1d_double(int n,
                                double *x,
                                double *y,
                                int nn,
                                double *xx,
                                double *yy,
                                int deriv);

    void deriv_pchip1d_double(int n,
                              double *x,
                              double *y,
                              int nn,
                              double *xx,
                              double *yy,
                              int deriv);

    void pchip1d_double(int n,
                        double *x,
                        double *y,
                        int nn,
                        double *xx,
                        double *yy);

    void quintic1d_double(int n,
                          double *x,
                          double *y,
                          int nn,
                          double *xx,
                          double *yy);

}

using namespace std;

// 1D linear spline interpolation
void lspline1d(int n,
               float *x,
               float *y,
               int nn,
               float *xx,
               float *yy)
{

    // Using cubic spline
    SplinesLoad::LinearSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    for (int i = 0; i < nn; i++)
    {
        yy[i] = sp(xx[i]);
    }
}

// 1D cubic spline interpolation
void cspline1d(int n,
               float *x,
               float *y,
               int nn,
               float *xx,
               float *yy)
{

    // Using cubic spline
    SplinesLoad::CubicSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    for (int i = 0; i < nn; i++)
    {
        yy[i] = sp(xx[i]);
    }
}

void deriv_cspline1d(int n,
                     float *x,
                     float *y,
                     int nn,
                     float *xx,
                     float *yy,
                     int deriv)
{

    // Using cubic spline
    SplinesLoad::CubicSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    switch (deriv)
    {
        case (1):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.D(xx[i]);
            }
        case (2):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.DD(xx[i]);
            }
        case (3):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.DDD(xx[i]);
            }
        default:
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.D(xx[i]);
            }
    }
}

void deriv_pchip1d(int n,
                   float *x,
                   float *y,
                   int nn,
                   float *xx,
                   float *yy,
                   int deriv)
{

    // Using cubic spline
    SplinesLoad::PchipSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    switch (deriv)
    {
        case (1):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.D(xx[i]);
            }
        case (2):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.DD(xx[i]);
            }
        case (3):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.DDD(xx[i]);
            }
        default:
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.D(xx[i]);
            }
    }
}

// 1D PCHIP interpolation
void pchip1d(int n,
             float *x,
             float *y,
             int nn,
             float *xx,
             float *yy)
{

    // Using PCHIP spline
    SplinesLoad::PchipSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    for (int i = 0; i < nn; i++)
    {
        yy[i] = sp(xx[i]);
    }
}

// 1D quntic interpolation
void quintic1d(int n,
               float *x,
               float *y,
               int nn,
               float *xx,
               float *yy)
{

    // Using quintic spline
    SplinesLoad::QuinticSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    for (int i = 0; i < nn; i++)
    {
        yy[i] = sp(xx[i]);
    }
}


// 1D linear spline interpolation
void lspline1d_double(int n,
                      double *x,
                      double *y,
                      int nn,
                      double *xx,
                      double *yy)
{

    // Using cubic spline
    SplinesLoad::LinearSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    for (int i = 0; i < nn; i++)
    {
        yy[i] = sp(xx[i]);
    }
}

// 1D cubic spline interpolation
void cspline1d_double(int n,
                      double *x,
                      double *y,
                      int nn,
                      double *xx,
                      double *yy)
{

    // Using cubic spline
    SplinesLoad::CubicSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    for (int i = 0; i < nn; i++)
    {
        yy[i] = sp(xx[i]);
    }
}

void deriv_cspline1d_double(int n,
                            double *x,
                            double *y,
                            int nn,
                            double *xx,
                            double *yy,
                            int deriv)
{

    // Using cubic spline
    SplinesLoad::CubicSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    switch (deriv)
    {
        case (1):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.D(xx[i]);
            }
        case (2):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.DD(xx[i]);
            }
        case (3):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.DDD(xx[i]);
            }
        default:
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.D(xx[i]);
            }
    }
}

// 1D PCHIP interpolation
void pchip1d_double(int n,
                    double *x,
                    double *y,
                    int nn,
                    double *xx,
                    double *yy)
{

    // Using PCHIP spline
    SplinesLoad::PchipSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    for (int i = 0; i < nn; i++)
    {
        yy[i] = sp(xx[i]);
    }
}

void deriv_pchip1d_double(int n,
                          double *x,
                          double *y,
                          int nn,
                          double *xx,
                          double *yy,
                          int deriv)
{

    // Using cubic spline
    SplinesLoad::PchipSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    switch (deriv)
    {
        case (1):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.D(xx[i]);
            }
        case (2):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.DD(xx[i]);
            }
        case (3):
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.DDD(xx[i]);
            }
        default:
            for (int i = 0; i < nn; i++)
            {
                yy[i] = sp.D(xx[i]);
            }
    }
}

// 1D quntic interpolation
void quintic1d_double(int n,
                      double *x,
                      double *y,
                      int nn,
                      double *xx,
                      double *yy)
{

    // Using quintic spline
    SplinesLoad::QuinticSpline sp;

    // Original data
    for (int i = 0; i < n; i++)
    {
        sp.pushBack(x[i], y[i]);
    }

    // Build the interpolation spline
    sp.build();

    // Compute spline values at new locations
    for (int i = 0; i < nn; i++)
    {
        yy[i] = sp(xx[i]);
    }
}
