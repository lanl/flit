# FLIT: A Generic Fortran Library based on Interfaces and Templates

## 1 Summary

We develop a generic Fortran library, `FLIT` (**F**ortran **L**ibrary based on **I**nterfaces and **T**emplates), to provide a number of useful functionalities for computational geophysics and beyond. These functionalities include single- and multidimensional array manipulation functions and subroutines, flexible parameter reading from text files or command-line arguments, signal and image filtering and processing, integral transforms, interpolation, statistical functions, and so on. These functionalities are not intrinsically available in the current Fortran standard. The most notable feature of `FLIT` is that it provides user-friendly interfaces for similar functionalities with different data types, demonstrating the possibility of accelerating scientific application development with modern Fortran. The package can serve as a robust Fortran library for developing sophisticated scientific applications in computational physics, computational and applied geophysics, signal and image processing, and so on.

## 2 Statement of need

Since its birth in 1957, Fortran has been a programming language of choice for developing scientific and engineering applications because of its simplicity, its efficiency, and, in particular, its capability to handle multidimensional arrays conveniently. Although newer programming languages such as Python, Julia, and modern C++ are gradually superseding it, modern Fortran (Fortran 90 and later) still anchors a wide variety of computationally intensive applications. However, in many cases, one still needs to develop low-level functions or subroutines from scratch to perform basic array operations. Many efforts have emerged in the past decade to provide a so-called standard library for Fortran (e.g., [stdlib](https://github.com/fortran-lang/stdlib)), covering array manipulation, string functions, hash tables, and many others, with varying levels of success and community acceptance.

In computational geophysics and many similar fields, one needs to routinely manipulate various types of arrays to achieve efficient manipulation, interpolation, filtering, and other processing of digitized signals and images. Unlike users of Python or MATLAB, many Fortran users still struggle to find a consistent, flexible, and efficient library that provides these functionalities. One of the most critical issues of Fortran is its lack of convenient generic programming (e.g., templates in C++) to accommodate different types of input variables. There are several choices for generic programming in Fortran (e.g., [fypp](https://github.com/aradi/fypp)), with different degrees of maturity and ease of use.

`FLIT` builds a generic Fortran library for scientific applications, especially computational physics and geophysics, by exploiting interfaces and templates. It contains modules for array allocation and manipulation, signal and image interpolation and filtering, flexible parameter reading, statistical functions, Fourier transforms, linear algebra, and so on. In particular, we design user-friendly interfaces for these functions and subroutines, so that a user can easily build sophisticated computational functionalities based on the modules. Unlike `fypp`, `FLIT` uses a more traditional approach based on preprocessor macros (`#define` and `#undef`) and file inclusion (`#include`) to achieve generic programming for various types of functionality. We emphasize that `FLIT` does not attempt to replace the current community effort to construct a Fortran standard library or anything alike. Nor does it attempt to be a *library of everything*, which is probably too ambitious for modern scientific research with highly diverse needs from different fields.

`FLIT` is designed for researchers in computational geophysics, and in a much broader scope, who need to routinely manipulate single- and multidimensional arrays, use signal and image processing, interpolation, complex parameter input, and so on. We have used it as a base library to develop sophisticated applications for seismic and electromagnetic wave modeling, traveltime tomography, full-wavefield imaging, and full-waveform inversion, as well as to generate high-fidelity synthetic data and images for training machine learning models. This work has led to a number of publications, including Gao et al. (2023), Gao (2024), and Gao and Modrak (2024). The generic functionalities of `FLIT` will likely enable more scientific applications, not only in computational geophysics but also in other computational fields.

## 3 Design and methodology

While `FLIT` contains many modules that cover a variety of needs of Fortran users, we highlight the following modules, which may be particularly relevant to computational geophysics applications. Wherever possible, we implement the functionalities in `FLIT` with OpenMP (for shared-memory parallelism) and/or the Message Passing Interface (MPI, for distributed-memory parallelism) to improve computational efficiency and scalability. [Section 4](#4-functions-in-flit) briefly lists the functions of all modules.

### 3.1 Array manipulation

We design a number of basic array allocation and manipulation functions and subroutines, and we provide a unified interface for each of them across different data types, extending the intrinsic procedures of the current Fortran standard.

- Array allocation: The `allocate` statement in the current Fortran standard cannot be applied to an already allocated array before the array is deallocated. This restriction can make code lengthy. `FLIT` provides a number of array allocation procedures for safe deallocation and allocation, e.g.,
  - `alloc_array(w, [lbound1, ubound1, ...], pad=1, source=s)` – a subroutine that combines `deallocate` and `allocate` based on specific lower and upper bounds, and also pads the array if needed. When `w` is already allocated, `alloc_array` first deallocates the array and then allocates it, ensuring memory safety;
  - `zeros(n1, n2, ...)` – a function that allocates an array of zeros. Through Fortran's automatic reallocation on assignment, the resulting array adjusts to the new dimensions if it is already allocated;
  - `ones(n1, n2, ...)` – a function that allocates an array of ones;
  - `const(n1, n2, ..., val)` – a function that allocates an array filled with the constant `val`;
  - `trues(n1, n2, ...)`, `falses(n1, n2, ...)` – functions that allocate an array of Fortran `.true.` or `.false.`;
  - `zeros_like`, `ones_like`, `const_like`, `trues_like`, `falses_like` – functions that allocate an array like another array (i.e., without specifying the dimensions);
  - `random(n1, n2, ..., dist, seed)` – a function that allocates an array filled with random values, where `dist` is the random distribution (`uniform`, `normal`, `cauchy`, `poisson`, or `exp`), and `seed` is the seed for reproducibility. Without `seed`, the random values differ every time the array is initialized, even in OpenMP or MPI processes that run simultaneously.
- Regularly spaced arrays: `FLIT` also provides convenient interfaces for generating regularly spaced arrays:
  - `linspace(a, b, n)` – generates a regularly spaced array of $n$ points from $a$ to $b$, i.e.,

    $$x_i = a + (i-1)\frac{b-a}{n-1}, \quad i = 1, 2, \ldots, n,$$

    where $x_i$ is the $i$th element of the resulting array. The resulting array has $n$ elements, with its first and last elements equal to $a$ and $b$, respectively;
  - `regspace(a, d, b)` – generates a regularly spaced array from $a$ with an interval of $d$ up to $b$, i.e.,

    $$x_i = a + (i-1)d, \quad i = 1, 2, \ldots, \quad x_i \le b.$$

    The first element of the resulting array is $a$, the interval between two consecutive elements is $d$, and the last element is the largest $x_i$ that does not exceed $b$. The number of elements is available through `size(x)`.

  Similar functions include `logspace`, which generates a regularly spaced array on a logarithmic scale, and `meshgrid`, which returns the coordinates of the points of a regular grid in one to three dimensions.
- Array manipulation: `FLIT` provides a number of array manipulation procedures, e.g.,
  - `crop` – crops an array;
  - `flatten` – flattens a multidimensional array to a 1D array;
  - `flip` – flips an array along one or more axes;
  - `pad` – pads an array at its edges;
  - `permute` – permutes the dimensions of an array, with the target order of dimensions specified as an integer, e.g., `order=321` and `order=2143`. `FLIT` interprets the integer `order` digit by digit as the target order. For instance, `321` means that after permutation, the third dimension moves to the first dimension, the second dimension stays the same, and the first dimension moves to the third dimension;
  - `rot90` – rotates an array by 90 degrees;
  - `tile` – repeats an array along one or more axes;
  - `rescale` – linearly rescales an array to a specific range;
  - `taper` – tapers an array with a window function. `FLIT` implements the step, linear, Parzen, Welch, sine, power-of-sine, Hann, Hamming, Blackman, Nuttall, Blackman–Nuttall, Blackman–Harris, Kaiser, and Gauss windows. [Wikipedia](https://en.wikipedia.org/wiki/Window_function) gives the definitions of these window functions.

### 3.2 Flexible parameter reading

The current Fortran standard provides a namelist-based parameter input method, which was a common vendor extension before it became standard in Fortran 90. However, in addition to requiring a predefined namelist in the code, this method falls short when a user needs more sophisticated parameter input, e.g., reading an unknown number of values into an array parameter, or reading a value into a parameter under certain conditions.

`FLIT` reads parameters from three types of sources: an ASCII (American Standard Code for Information Interchange) text file, command-line arguments, and an array of strings. Currently, `FLIT` reads integers (kinds 2, 4, and 8), single- and double-precision real numbers, single- and double-precision complex numbers, strings, and logical values. The corresponding suffixes `#type` (used below) are `int` (for all three integer kinds), `float`, `double`, `complex`, `dcomplex`, `string`, and `logical`.

For parameter input from an ASCII file, `FLIT` provides the following interfaces:

- Read a single value into a variable:

  ```fortran
  call readpar_#type(parameter_file, parameter_name, target_variable, default_value, required)
  ```

  where `parameter_file` is the ASCII file that holds the parameters, `parameter_name` is the name of the parameter, `target_variable` is the variable that stores the parameter value, `default_value` is the value of the parameter if it is not in the parameter file or has an empty value (e.g., `n =`), and `required` is an optional logical argument that states whether the parameter is required. If a required parameter is not in the parameter file, the program reports an error and stops.
- Read multiple values into an allocatable array:

  ```fortran
  call readpar_n#type(parameter_file, parameter_name, target_variable, default_value, required, separator)
  ```

  The parameter file can give the values as `var = 1, 2, 3, 4` for an integer or float parameter, `var = abcd, efgh, x, y, z` for a string parameter, or `var = y, n, n, y` for a logical parameter. The optional argument `separator` sets the delimiter between the values (the default is a comma); for example, `separator=';'` reads `var = 1; 2; 3; 4`.
- Read a value that depends on a condition:

  ```fortran
  call readpar_x#type(parameter_file, parameter_name, target_variable, default_value, condition, required)
  ```

  For instance, suppose that the parameter file specifies `par = 0:0, 10:1, 11~20:2`, which means that `par = 0` when `condition = 0`, `par = 1` when `condition = 10`, and `par = 2` when `condition` is between 11 and 20, where `~` denotes a range. Then `call readpar_xfloat(parameter_file, 'par', par, 1.0, 3.0)`, i.e., with `condition = 3.0`, returns `par = 0.3`, because `FLIT` linearly interpolates between the known values 0 and 1 to obtain `par` at `condition = 3`, which the parameter file does not specify. For integer, logical, and string parameters, where linear interpolation does not apply, `FLIT` uses nearest-neighbor interpolation instead. For instance, with the same specification, `call readpar_xint(parameter_file, 'par', par, 1, 3.0)` returns `par = 0`, because 3.0 is closer to 0 than to 10. Halliburton's Landmark ProMAX (see, for example, the [ProMAX documentation](https://esd.halliburton.com/support/LSM/GGT/ProMAXSuite/ProMAX/5000/5000_8/Help/promax/user_int.pdf)), a well-known seismic data processing package that newer Halliburton products have superseded, inspired this parameter input method. In seismic processing, one often needs to specify the value of a parameter based on the value of another variable (i.e., the condition here), and sometimes the other variable is given as a range, similar to the `~` notation here. We believe that this parameter reading method is a useful feature for developing computational applications, yet it is not available in the current Fortran standard.

Similarly, `FLIT` provides `getpar_#type`, `getpar_n#type`, and `getpar_x#type` to read parameters from command-line arguments, and `parsepar_#type`, `parsepar_n#type`, and `parsepar_x#type` to read parameters from an array of strings. Their syntax is similar to that of `readpar`, except for the source: command-line arguments for `getpar`, and an array of strings for `parsepar`. On the command line, write each parameter as `name=value` without spaces around `=`.

It is worth mentioning that `FLIT` reads complex numbers flexibly. For instance, `FLIT` reads both `var = 1.0 + 2.0i` and `var = 2i + 1` as the complex number $1 + 2i$, where $i$ is the imaginary unit. This notation is more concise and natural than Fortran's native `(1.0, 2.0)` syntax.

### 3.3 Image and signal processing

The image and signal filters in `FLIT` include the following.

- Gaussian smoothing: Gaussian filtering is one of the most frequently used filters in image processing and in many geophysical applications. `FLIT` implements four methods of Gaussian filtering: the Deriche recursive filter (Deriche, 1990), filtering in the discrete cosine transform (DCT) domain, the recursive spectral filter (Sugimoto and Kamata, 2013), and kernel convolution based on the Fourier transform;
- Laplacian filtering: Applies the discrete Laplacian operator, which suppresses low-wavenumber or low-frequency components. [Wikipedia](https://en.wikipedia.org/wiki/Discrete_Laplace_operator) gives the details;
- Mean filtering: Smooths an array with the mean of a moving window;
- Median filtering: Smooths an array with the median of a moving window;
- LOWESS filtering: Locally weighted scatterplot smoothing (LOWESS) is a popular nonparametric tool in regression analysis, and it can also act as a filter that smooths out outliers in data with local polynomials. We implement LOWESS filters for 1D, 2D, and 3D data. [Wikipedia](https://en.wikipedia.org/wiki/Local_regression) gives the details;
- Fourier filtering: Filters a signal $s(t)$ by

  $$s'(t) = \Re\left[\mathcal{F}^{-1}\left[K(\omega)\mathcal{F}[s(t)]\right]\right],$$

  where $t$ is time (or space), $\omega$ is the angular frequency (or wavenumber), $s'(t)$ is the filtered signal, $\mathcal{F}$ and $\mathcal{F}^{-1}$ are the forward and inverse Fourier transforms, $\Re$ takes the real part, and $K(\omega)$ is a real, zero-phase amplitude response that the user specifies at a set of frequencies (or wavenumbers). Different choices of $K(\omega)$ result in lowpass, highpass, bandpass, bandstop, and other filters. In 2D and 3D, the filter is the product of the 1D amplitude responses along the axes;
- Dip filtering: Filters a multidimensional image based on its local dip. `FLIT` computes the local dip, or slope, from the eigenvectors of the Gaussian-smoothed structure tensor (Bakker, 2002);
- Windowed RMS filtering: Balances the amplitude of a single- or multidimensional array with the windowed root mean square (RMS). In 1D,

  $$s'_i = \frac{s_i}{\sqrt{\frac{1}{2w+1}\left(\sum_{j=-w}^{w} s_{i+j}^2 + \epsilon\right)}},$$

  where $s_i$ and $s'_i$ are the $i$th samples of the input and output arrays, $w$ is the half length of the window, and $\epsilon$ is a small positive number that prevents division by zero (by default, 0.01 times the largest $s_{i+j}^2$ in the window);
- Nonlinear anisotropic diffusion filtering (AnDF): We implement a fast explicit diffusion (FED)-based AnDF solver (Weickert et al., 2016) for both 2D and 3D cases. The solver relies on structure tensors and includes some improvements over the original algorithm. We also provide MPI versions of the 2D and 3D solvers to accelerate the computation for large images;
- Total variation and total generalized $p$-variation filtering: The total variation (TV) filter is an important filter in image processing and, more recently, in seismic tomography and full-waveform inversion. We implement isotropic TV filtering and total generalized $p$-variation (TGpV) filtering, where $p$ is the exponent of the $\ell_p$ norm. For isotropic TV, we follow the algorithm of Goldstein and Osher (2009). For TGpV, we implement the algorithms of Knoll et al. (2011) and Gao and Huang (2019). `FLIT` also provides MPI versions of the 2D and 3D TGpV filters.

In addition, we implement single- and multidimensional convolutions with both a direct (discrete) algorithm, `convd`, and a Fourier-transform-based algorithm, `conv`. The latter has a unified interface for 1D, 2D, and 3D arrays,

$$s' = \mathrm{conv}(s, K, \mathrm{method}),$$

where $s$ is the input array, $K$ is the discrete kernel specified by the user, and $s'$ is the resulting array. When `method = 'same'`, $s'$ has the same size as $s$.

### 3.4 Miscellaneous functionalities

`FLIT` also contains many functionalities that do not fit in the above three categories, including statistical functions, integral transforms (Fourier, Hilbert, DCT, etc.), differentiation operators, and geometry (points, lines, surfaces, rotation, etc.). None of these functionalities is part of the current Fortran standard. [Section 4](#4-functions-in-flit) lists them, and the source code gives the details.

## 4 Functions in FLIT

The tables below briefly list the public functions, subroutines, and types of `FLIT`, grouped by area. A program can access all of them with `use libflit`, which also includes the k-d tree module `kdtree_mod` and, when `FLIT` is built with HDF5, the HDF5 module `fh5`. Most names are generic interfaces that accept several data types and array dimensions.

### 4.1 Arrays

| Module | Purpose | Functions, subroutines, and types |
|---|---|---|
| `libflit_array` | Array allocation and regularly spaced arrays | `alloc_array`, `zeros`, `ones`, `const`, `trues`, `falses`, `zeros_like`, `ones_like`, `const_like`, `trues_like`, `falses_like`, `eye`, `diag`, `linspace`, `regspace`, `logspace` |
| `libflit_array_operation` | Array manipulation | `crop`, `crop_array`, `flatten`, `flip`, `pad`, `pad_array`, `permute`, `permute_array`, `rot90`, `tile`, `rescale`, `slice`, `mask`, `ifelse`, `binarize`, `adjust`, `as_scalar`, `cross2`, `cross3`, `ndigits`, `breakint`, `any_in`, `all_in`, `remove_any_in` |
| `libflit_array_extension` | Named arrays, each stored with a name, an ID, and its dimensions | the types `meta_array1_#type` to `meta_array4_#type`, where `#type` is `int`, `real`, `double`, `complex`, `dcomplex`, or `logical`, and `get_meta_array`, `get_meta_array_core`, `set_meta_array_core` |
| `libflit_taper` | Window functions and tapering | `taper`, `taper_window`, `window_function` |
| `libflit_sort` | Sorting | `sort`, `sort_index`, `mergesort_index` |
| `libflit_unique` | Unique elements | `unique` |

### 4.2 Parameter reading

| Module | Purpose | Functions, subroutines, and types |
|---|---|---|
| `libflit_readpar` | Parameter reading from files, command-line arguments, or arrays of strings | `readpar_#type`, `readpar_n#type`, `readpar_x#type`, `getpar_#type`, `getpar_n#type`, `getpar_x#type`, `parsepar_#type`, `parsepar_n#type`, `parsepar_x#type`, `checkpar`, `checkarg` |

Here `#type` is `int`, `float`, `double`, `complex`, `dcomplex`, `logical`, or `string`. The subroutines `checkpar` and `checkarg` stop the program if a parameter file or the command line contains a parameter name that is not in a given list of valid names.

### 4.3 Filters

| Module | Purpose | Functions, subroutines, and types |
|---|---|---|
| `libflit_gaussfilt` | Gaussian filter | `gauss_filt` |
| `libflit_laplacefilt` | Laplacian filter | `laplace_filt` |
| `libflit_meanfilt` | Mean filter | `mean_filt` |
| `libflit_medianfilt` | Median filter | `median_filt` |
| `libflit_lowessfilt` | LOWESS filter | `lowess_filt` |
| `libflit_fourierfilt` | Fourier filters | `fourier_filt`, `fourier_filter`, `fourier_smooth`, `fourier_sharpen` |
| `libflit_dip` | Dip filter | `local_dip`, `dip_filt` |
| `libflit_balancefilt` | Windowed RMS filter | `balance_filt` |
| `libflit_andffilt` | Nonlinear anisotropic diffusion filter | `andf_filt`, `andf_filt_mpi`, `compute_structure_tensor`, `compute_structure_tensor_eigens`, and the type `andf_param` |
| `libflit_tvfilt` | TV and TGpV filters | `tv_filt`, `sparse_tv_filt`, `tgpv_filt`, `tgpv_filt_mpi`, `soft_shrinkage` |
| `libflit_iirfilt` | Infinite impulse response (IIR) filters | `iir_filt` |

### 4.4 Transforms, interpolation, and calculus

| Module | Purpose | Functions, subroutines, and types |
|---|---|---|
| `libflit_transform` | Integral transforms and convolution | `fft`, `ifft`, `fftd`, `ifftd`, `fourier_transform`, `inverse_fourier_transform`, `fftshift`, `fft_omega`, `fft_deriv`, `cosine_transform`, `inverse_cosine_transform`, `dct`, `idct`, `hilbert`, `hilbert_transform`, `envelope`, `envelope_transform`, `phase_shift`, `time_shift`, `conv`, `convd`, `deconv`, `next_power_2`, `next_power_235`, `next_power_2357`, `next_power_235711` |
| `libflit_spectrum` | Time-frequency analysis | `gst`, `igst`, `gst_simple`, `igst_simple`, `gabor`, `igabor`, `gabor_simple`, `igabor2`, `tfspec`, `instant_phase` |
| `libflit_interp` | Interpolation | `interp`, `interp_to`, `interp_like`, `resample`, `ginterp`, `meshgrid`, `inpaint`, `point_interp_linear`, `point_interp_barycentric`, and the logical switch `mba_verbose` |
| `libflit_fit` | Polynomial fitting | the types `polyfit1` and `polyfit2` |
| `libflit_calculus` | Derivatives and integrals | `deriv`, `integ`, `cumsum`, `taylor_fd_coefs` |
| `libflit_specialfunc` | Special functions | `sinc`, `dsinc`, `dirac_delta`, `bessel_i0`, `bessel_k0`, `factorial`, `binomial`, `ei`, `expei`, `eone`, `cgamma` |

### 4.5 Statistics, random numbers, and linear algebra

| Module | Purpose | Functions, subroutines, and types |
|---|---|---|
| `libflit_statistics` | Statistics and correlation | `mean`, `median`, `var`, `std`, `covar`, `covar_matrix`, `histogram`, `plot_histogram`, `xcorr`, `acorr`, `xcorrd`, `acorrd`, `lxcorr`, `xcorr_coef`, `gaussian`, `cauchy`, `kernel_triangular`, `kernel_bisquare`, `kernel_trisquare`, `kernel_tricube`, `kernel_epanechnikov`, `kernel_cosine` |
| `libflit_random`, `libflit_random_extra` | Random numbers | `rand`, `drand`, `irand`, `random`, `drandom`, `irandom`, `random_order`, `random_permute`, `random_mask`, `random_string`, `random_pdf` |
| `libflit_linear_algebra` | Linear algebra | `det`, `inv`, `trace`, `solve`, `lsqsolve`, `solve_band`, `svd`, `eigen`, `eigen_symm2x2`, `eigen_symm3x3`, `spectral_radius`, `matx`, `diagx`, `xdiag`, `toeplitz_matrix`, `hankel_matrix`, `vandermonde_matrix` |

### 4.6 Geometry, clustering, and data structures

| Module | Purpose | Functions, subroutines, and types |
|---|---|---|
| `libflit_geometry` | Computational geometry | `convexhull`, `polygon_contains_point`, `rotate_point`, `rotation_matrix`, `point_distance_to_line`, `point_signed_distance_to_line`, `point_distance_to_line_segment`, `distance_between_line_segments`, `gaussian_curvature`, `spherical_to_cartesian`, `cartesian_to_spherical`, `fit_curve`, `fit_surface` |
| `libflit_clustering` | Clustering | `hdbscan_clustering` |
| `libflit_heap` | Heaps | the types `heap_int`, `heap_float`, and `heap_double`, with `init`, `insert`, `peek`, `pop`, `reheap`, `size`, and `free` |
| `kdtree_mod` | k-d tree | the type `kdtree`, with `build`, `query_nn`, `query_knn`, `query_ball`, `num_points`, `num_dims`, `built`, and `free` |

### 4.7 Strings, files, and utilities

| Module | Purpose | Functions, subroutines, and types |
|---|---|---|
| `libflit_string` | String manipulation and conversion | `num2str`, `array_to_string`, `tidy`, `to_upper`, `to_lower`, `split_string`, `substring`, `reverse_string`, `center_substring`, `count_substring`, `count_words`, `get_string_after`, `remove_string_after`, and the `extract_*` functions that convert strings to values |
| `libflit_io` | Binary array input and output | `input_array`, `output_array`, `stdin_array`, `stdout_array`, `load` |
| `libflit_filedir` | Files and directories | `file_exists`, `directory_exists`, `is_directory`, `get_file_size`, `get_stdin_size`, `copy_file`, `move_file`, `delete_file`, `list_files`, `count_files`, `count_lines`, `count_nonempty_lines`, `make_directory`, `move_directory`, `delete_directory`, `get_current_directory`, `get_real_path`, `get_file_directory`, `get_basename`, `get_extension`, `strip_extension` |
| `libflit_date_time` | Date and time | `date_time_string`, `date_time_stamp`, `date_time_compact`, `print_date_time`, `time_duration` |
| `libflit_utility` | General utilities | `warn`, `clip`, `swap`, `round`, `nice`, `rov`, `order_of_magnitude`, `get_integer_digits`, `index_first_nonzero`, `equivalent`, `return_normal`, `print_array`, `get_hostname` |
| `libflit_error` | Assertions | `assert`, `remind` |
| `libflit_constants` | Constants | named constants such as `const_pi`, `const_deg2rad`, `const_lightc`, `float_huge`, `double_tiny`, `float_nan`, and `char_degree`, and the function `nan` |

### 4.8 Parallel computing and HDF5

| Module | Purpose | Functions, subroutines, and types |
|---|---|---|
| `libflit_mpicomm` | MPI communication | `mpistart`, `mpiend`, `mpibarrier`, `mpistop`, `bcast`, `bcast_array`, `reduce`, `reduce_array`, `allreduce`, `allreduce_array`, `commute_array`, `domain_decomp_regular`, `global_min`, `global_max`, `global_and`, `global_or` |
| `libflit_mpicomm_group` | MPI communication within a group of processes | the routines of `libflit_mpicomm` with the suffix `_group` (e.g., `bcast_group` and `commute_array_group`), and `group_min`, `group_max`, `group_and`, `group_or` |
| `libflit_domain_decomposition` | Domain decomposition | `divide_domain`, `cut`, `split_list` |
| `fh5` | HDF5 input and output | `fh5_open`, `fh5_close`, `fh5_finalize`, `fh5_write`, `fh5_read`, `fh5_write_str`, `fh5_read_str`, `fh5_write_attr`, `fh5_read_attr`, `fh5_print_attr`, `fh5_exists`, `fh5_ndim`, `fh5_shape`, `fh5_dtype`, `fh5_list`, `fh5_print_tree`, `fh5_write_chunk`, `fh5_delete`, `fh5_copy`, `fh5_fid` |

## 5 Acknowledgments

This work was supported by the Los Alamos National Laboratory (LANL) Laboratory Directed Research and Development (LDRD) program. LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration (NNSA) of the U.S. Department of Energy (DOE) under Contract No. 89233218CNA000001. This research used high-performance computing resources provided by LANL's Institutional Computing program.

## References

Bakker, P., 2002, Image structure analysis for seismic interpretation: Ph.D. dissertation, Delft University of Technology, [citeseerx.ist.psu.edu](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=ab677525a49b902a42e3e3b42e779dcbcb957b1f).

Deriche, R., 1990, Fast algorithms for low-level vision: IEEE Transactions on Pattern Analysis and Machine Intelligence, **12**, no. 1, 78–87, doi: [10.1109/34.41386](https://doi.org/10.1109/34.41386).

Gao, K., 2024, Iterative multitask learning and inference from seismic images: Geophysical Journal International, **236**, no. 1, 565–592, doi: [10.1093/gji/ggad424](https://doi.org/10.1093/gji/ggad424).

Gao, K., C. Donahue, B. G. Henderson, and R. T. Modrak, 2023, Deep-learning-guided high-resolution subsurface reflectivity imaging with application to ground-penetrating radar data: Geophysical Journal International, **233**, no. 1, 448–471, doi: [10.1093/gji/ggac468](https://doi.org/10.1093/gji/ggac468).

Gao, K., and L. Huang, 2019, Acoustic- and elastic-waveform inversion with total generalized p-variation regularization: Geophysical Journal International, **218**, no. 2, 933–957, doi: [10.1093/gji/ggz203](https://doi.org/10.1093/gji/ggz203).

Gao, K., and R. T. Modrak, 2024, Machine learning inference of random medium properties: IEEE Transactions on Geoscience and Remote Sensing, **62**, 1–13, doi: [10.1109/TGRS.2024.3367541](https://doi.org/10.1109/TGRS.2024.3367541).

Goldstein, T., and S. Osher, 2009, The split Bregman method for L1-regularized problems: SIAM Journal on Imaging Sciences, **2**, no. 2, 323–343, doi: [10.1137/080725891](https://doi.org/10.1137/080725891).

Knoll, F., K. Bredies, T. Pock, and R. Stollberger, 2011, Second order total generalized variation (TGV) for MRI: Magnetic Resonance in Medicine, **65**, no. 2, 480–491, doi: [10.1002/mrm.22595](https://doi.org/10.1002/mrm.22595).

Sugimoto, K., and S.-i. Kamata, 2013, Fast Gaussian filter with second-order shift property of DCT-5: 2013 IEEE International Conference on Image Processing, 514–518, doi: [10.1109/ICIP.2013.6738106](https://doi.org/10.1109/ICIP.2013.6738106).

Weickert, J., S. Grewenig, C. Schroers, and A. Bruhn, 2016, Cyclic schemes for PDE-based image analysis: International Journal of Computer Vision, **118**, no. 3, 275–299, doi: [10.1007/s11263-015-0874-1](https://doi.org/10.1007/s11263-015-0874-1).
