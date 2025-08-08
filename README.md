# Description
**FLIT: A Generic Fortran Library based on Interfaces and Templates**

We develop a generic Fortran library, `FLIT` (Fortran Library based on Interfaces and Templates), to provide a number of useful functionalities for computational geophysics and beyond. These functionalities include several single/multi-dimensional array manipulation functions/subroutines, flexible parameter reading from textual file or command line arguments, signal and image filtering and processing, integral transforms, interpolation, statistical functions, and so on. These functionalities are not intrinsically available in the current Fortran standard. The most notable feature of `FLIT` is that we provide user-friendly interfaces for similar functionalities with different data types, attempting to demonstrate the possibility of accelerating scientific application development using modern Fortran. The package can be used as a robust Fortran library for developing sophisticated scientific applications in the fields of computational physics, computational and applied geophysics, signal and image processing, and so on. 

The work is supported by Los Alamos National Laboratory (LANL) Laboratory Directory Research and Development (LDRD) project 20240322ER. LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration (NNSA) of the U.S. Department of Energy (DOE) under Contract No. 89233218CNA000001. The research used high-performance computing resources provided by LANL's Institutional Computing program. 

The work is under LANL open source approval reference O4767.

# Reference
Please refer to the document [LA-UR-24-26315](doc/doc_libflit.pdf) for details. 

# Requirement
Currently, `FLIT` can only be compiled with Intel's oneAPI HPC Toolkit 2024.2/2025.2 (including ifx, icx, icpx, and mpiifx compilers) on Linux platforms. You need to install [Intel's oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html) and set the paths using the script [set_intel_paths.sh](set_intel_paths.sh) included this repository. `FLIT` cannot be compiled with GNU's compiler suite (gcc and gfortran). 

# Use
The codes are written in Fortran, with some functions written in C++ interfacing to third-party libraries in [third_party](third_party). 

To install `FLIT`, 

```
cd src
make
```

The compiled `FLIT` will be at `lib` including module/submodule files and a static library file, `libflit.a`. 

To remake, 

```
cd src
make clean
make
```

We include several simple examples to use `FLIT` in `test`. To run the tests, 

```
cd test
bash test.sh
```

The [Makefile](src/Makefile) in the [test](test) directory can serve as an example on how to use `FLIT` in your code, including path inclusion and link of the compiled library. 

Third-party codes are included in [third_party](third_party) for the purpose of completeness. Please read [third_party/README](third_party/README) for details. 

# License
&copy; 2024 - 2025. Triad National Security, LLC. All rights reserved. 

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Author
Kai Gao, <kaigao@lanl.gov>

We welcome feedbacks, bug reports, and improvement ideas on `FLIT`. 

If you use this package in your research and find it useful, please cite it as

* Kai Gao, Ting Chen, 2024, FLIT: A Generic Fortran Library based on Interfaces and Templates, url: [github.com/lanl/flit](https://github.com/lanl/flit)
