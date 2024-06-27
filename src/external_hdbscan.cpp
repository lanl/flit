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

#include <iostream>
#include <stdio.h>
#include "hdbscan.hpp"

using namespace std;

// external interface
extern "C" {

    void hdbscan(int n,
                 int nd,
                 float *x,
                 int min_sample,
                 int min_cluster_size,
                 int metric,
                 float p,
                 int &ncluster,
                 int &nnoisy,
                 int *labels,
                 float *probs);
}


void hdbscan(int n,
             int nd,
             float *x,
             int min_sample,
             int min_cluster_size,
             int metric,
             float p,
             int &ncluster,
             int &nnoisy,
             int *labels,
             float *probs)
{

    Hdbscan hc;

    // The input data are organized as a 1D array tuple by tuple, i.e.,
    // (x1, y1, ..), (x2, y2, ..), ..., (xn, y_n, ...)
    for (int i = 0; i < n; i++)
    {
        vector<double> element;
        for (int j = 0; j < nd; j++)
        {
            element.push_back(double(x[i * nd + j]));
        }
        hc.dataset.push_back(element);
    }

    switch (metric)
    {
        case 1:
            hc.execute(min_sample, min_cluster_size, "Euclidean");
        case 2:
            hc.execute(min_sample, min_cluster_size, "Manhattan");
        case 3:
            hc.execute(min_sample, min_cluster_size, "Minkowski", p);
    }

    // Copy results out
    ncluster = hc.numClusters_;
    nnoisy = hc.noisyPoints_;

    for (int i = 0; i < n; i++)
    {
        labels[i] = hc.labels_[i];
        probs[i] = hc.membershipProbabilities_[i];
    }

}
