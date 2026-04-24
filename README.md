# Description
**FLIT: A Generic Fortran Library based on Interfaces and Templates**

`FLIT` (Fortran Library based on Interfaces and Templates) is a generic Fortran library mainly for
 - Single/multi-dimensional array manipulation
 - Flexible parameter reading from textual file or command line arguments
 - Signal and image filtering and processing
 - Integral transforms
 - Interpolation
 - Statistical functions
 - etc. 
 These functionalities are not natively available in the current Fortran standard. 
 
In `FLIT`, we provide user-friendly interfaces for similar functionalities with different data types, attempting to demonstrate the possibility of accelerating scientific application development using modern Fortran. The package can be used as a robust Fortran library for developing sophisticated scientific applications in the fields of computational physics, computational and applied geophysics, signal and image processing, and so on. 

The work is supported by Los Alamos National Laboratory (LANL) Laboratory Directory Research and Development (LDRD) project 20240322ER. LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration (NNSA) of the U.S. Department of Energy (DOE) under Contract No. 89233218CNA000001. The research used high-performance computing resources provided by LANL's Institutional Computing program. 

The work is under LANL open source approval reference O4767.

# Requirement
 - Platform: Linux
 - Compiler: [Intel's oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html)

# Use
To set Intel's compiler paths: 
```
bash set_intel_paths.sh
```

To install `FLIT`:
```
cd src
make
```
The compiled `FLIT` will be at `lib` including module/submodule files and a single static library file, `libflit.a`. 

To remake:
```
cd src
make clean
make
```

To run examples:
```
cd test
bash test.sh
```

The [Makefile](src/Makefile) in the [test](test) directory can serve as an example on how to use `FLIT` in your code, including path inclusion and link of the compiled library. 

Third-party codes are included in [third_party](third_party). Please refer to [third_party/README](third_party/README.md) for details. 

# License
&copy; 2024-2026. Triad National Security, LLC. All rights reserved. 

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Documentation
Under development. Please refer to the document [LA-UR-24-26315](doc/doc_libflit.pdf) for details. 

# Author
Kai Gao, <kaigao@lanl.gov>

We welcome feedbacks, bug reports, and improvement ideas on `FLIT`. 

If you use this package in your research and find it useful, please cite it as
* Kai Gao, Ting Chen, 2024, FLIT: A Generic Fortran Library based on Interfaces and Templates, url: [github.com/lanl/flit](https://github.com/lanl/flit)
