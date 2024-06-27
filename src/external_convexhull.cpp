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

#define CONVHULL_3D_USE_SINGLE_PRECISION /* (optional) */
#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"
#include "convhull_2d.h"

#include <iostream>
#include <stdio.h>

using namespace std;

// external interface
extern "C" {

    void convhull_2d(int n,
                     float *x,
                     float *y,
                     int &nf,
                     int *f);

    void convhull_3d(int n,
                     float *x,
                     float *y,
                     float *z,
                     int &nf,
                     int *f);
}


void convhull_2d(int n,
                 float *x,
                 float *y,
                 int &nf,
                 int *f)
{
    std::vector<IndexedPoint> p;
    for (size_t i = 0; i < n; i++)
    {
        p.push_back({i, x[i], y[i]});
    }

    auto hull = convex_hull(p);

    // Copy results out
    nf = hull.size();
    for (size_t i = 0; i < nf; i++)
    {
        f[i] = hull[i].i;
    }

}


void convhull_3d(int n,
                 float *x,
                 float *y,
                 float *z,
                 int &nf,
                 int *f)
{

    ch_vertex* vertices;
    vertices = (ch_vertex*)malloc(n * sizeof(ch_vertex));
    for (int i = 0; i < n; i++)
    {
        vertices[i].x = x[i];
        vertices[i].y = y[i];
        vertices[i].z = z[i];
    }

    int nFaces;
    int *faceIndices = NULL;
    convhull_3d_build(vertices, n, &faceIndices, &nFaces);

    // Copy results out
    nf = nFaces;
    for (int i = 0; i < nf * 3; i++)
    {
        f[i] = faceIndices[i];
    }

}
