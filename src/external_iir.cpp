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

#include "iir.h"

// external interface
extern "C" {
    // Low-pass
    void butterworth_lowpass(const float *x,
                             float *xx,
                             const int n,
                             const float d,
                             const float freq,
                             const int order);
    void chebyshev1_lowpass(const float *x,
                            float *xx,
                            const int n,
                            const float d,
                            const float freq,
                            const int order,
                            const float pass_ripple_db);
    void chebyshev2_lowpass(const float *x,
                            float *xx,
                            const int n,
                            const float d,
                            const float freq,
                            const int order,
                            const float stop_ripple_db);
    // High-pass
    void butterworth_highpass(const float *x,
                              float *xx,
                              const int n,
                              const float d,
                              const float freq,
                              const int order);
    void chebyshev1_highpass(const float *x,
                             float *xx,
                             const int n,
                             const float d,
                             const float freq,
                             const int order,
                             const float pass_ripple_db);
    void chebyshev2_highpass(const float *x,
                             float *xx,
                             const int n,
                             const float d,
                             const float freq,
                             const int order,
                             const float stop_ripple_db);
    // Band-pass
    void butterworth_bandpass(const float *x,
                              float *xx,
                              const int n,
                              const float d,
                              const float freq1,
                              const float freq2,
                              const int order);
    void chebyshev1_bandpass(const float *x,
                             float *xx,
                             const int n,
                             const float d,
                             const float freq1,
                             const float freq2,
                             const int order,
                             const float pass_ripple_db);
    void chebyshev2_bandpass(const float *x,
                             float *xx,
                             const int n,
                             const float d,
                             const float freq1,
                             const float freq2,
                             const int order,
                             const float stop_ripple_db);
    // Band-stop
    void butterworth_bandstop(const float *x,
                              float *xx,
                              const int n,
                              const float d,
                              const float freq1,
                              const float freq2,
                              const int order);
    void chebyshev1_bandstop(const float *x,
                             float *xx,
                             const int n,
                             const float d,
                             const float freq1,
                             const float freq2,
                             const int order,
                             const float pass_ripple_db);
    void chebyshev2_bandstop(const float *x,
                             float *xx,
                             const int n,
                             const float d,
                             const float freq1,
                             const float freq2,
                             const int order,
                             const float stop_ripple_db);
}

using namespace std;

// 1D Butterworth low-pass filtering
void butterworth_lowpass(const float *x,
                         float *xx,
                         const int n,
                         const float d,
                         const float freq,
                         const int order)
{

    // Create the filter structure for n-th order
    Iir::Butterworth::LowPass<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, freq);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }

}

// 1D Chebyshev-I low-pass filtering
void chebyshev1_lowpass(const float *x,
                        float *xx,
                        const int n,
                        const float d,
                        const float freq,
                        const int order,
                        const float pass_ripple_db)
{

    // Create the filter structure for n-th order
    Iir::ChebyshevI::LowPass<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, freq, pass_ripple_db);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }
}


// 1D Chebyshev-II low-pass filtering
void chebyshev2_lowpass(const float *x,
                        float *xx,
                        const int n,
                        const float d,
                        const float freq,
                        const int order,
                        const float stop_ripple_db)
{

    // Create the filter structure for n-th order
    Iir::ChebyshevII::LowPass<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, freq, stop_ripple_db);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }
}


// 1D Butterworth high-pass filtering
void butterworth_highpass(const float *x,
                          float *xx,
                          const int n,
                          const float d,
                          const float freq,
                          const int order)
{

    // Create the filter structure for n-th order
    Iir::Butterworth::HighPass<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, freq);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }

}

// 1D Chebyshev-I high-pass filtering
void chebyshev1_highpass(const float *x,
                         float *xx,
                         const int n,
                         const float d,
                         const float freq,
                         const int order,
                         const float pass_ripple_db)
{

    // Create the filter structure for n-th order
    Iir::ChebyshevI::HighPass<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, freq, pass_ripple_db);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }
}


// 1D Chebyshev-II high-pass filtering
void chebyshev2_highpass(const float *x,
                         float *xx,
                         const int n,
                         const float d,
                         const float freq,
                         const int order,
                         const float stop_ripple_db)
{

    // Create the filter structure for n-th order
    Iir::ChebyshevII::HighPass<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, freq, stop_ripple_db);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }
}


// 1D Butterworth band-pass filtering
void butterworth_bandpass(const float *x,
                          float *xx,
                          const int n,
                          const float d,
                          const float freq1,
                          const float freq2,
                          const int order)
{

    // Create the filter structure for n-th order
    Iir::Butterworth::BandPass<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, 0.5 * (freq1 + freq2), freq2 - freq1);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }

}

// 1D Chebyshev-I bandpass filtering
void chebyshev1_bandpass(const float *x,
                         float *xx,
                         const int n,
                         const float d,
                         const float freq1,
                         const float freq2,
                         const int order,
                         const float pass_ripple_db)
{

    // Create the filter structure for n-th order
    Iir::ChebyshevI::BandPass<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, 0.5 * (freq1 + freq2), freq2 - freq1, pass_ripple_db);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }

}


// 1D Chebyshev-II bandpass filtering
void chebyshev2_bandpass(const float *x,
                         float *xx,
                         const int n,
                         const float d,
                         const float freq1,
                         const float freq2,
                         const int order,
                         const float stop_ripple_db)
{

    // Create the filter structure for n-th order
    Iir::ChebyshevII::BandPass<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, 0.5 * (freq1 + freq2), freq2 - freq1, stop_ripple_db);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }

}


// 1D Butterworth band-pass filtering
void butterworth_bandstop(const float *x,
                          float *xx,
                          const int n,
                          const float d,
                          const float freq1,
                          const float freq2,
                          const int order)
{

    // Create the filter structure for n-th order
    Iir::Butterworth::BandStop<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, 0.5 * (freq1 + freq2), freq2 - freq1);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }

}

// 1D Chebyshev-I bandstopfiltering
void chebyshev1_bandstop(const float *x,
                         float *xx,
                         const int n,
                         const float d,
                         const float freq1,
                         const float freq2,
                         const int order,
                         const float pass_ripple_db)
{

    // Create the filter structure for n-th order
    Iir::ChebyshevI::BandStop<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, 0.5 * (freq1 + freq2), freq2 - freq1, pass_ripple_db);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }

}


// 1D Chebyshev-II bandstopfiltering
void chebyshev2_bandstop(const float *x,
                         float *xx,
                         const int n,
                         const float d,
                         const float freq1,
                         const float freq2,
                         const int order,
                         const float stop_ripple_db)
{

    // Create the filter structure for n-th order
    Iir::ChebyshevII::BandStop<> f;

    // Compute the the coefficients
    f.setup(order, 1.0 / d, 0.5 * (freq1 + freq2), freq2 - freq1, stop_ripple_db);

    for (int i = 0; i < n; i++)
    {
        xx[i] = f.filter(x[i]);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        xx[i] = f.filter(xx[i]);
    }

}



