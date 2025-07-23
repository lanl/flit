
# Paths
bindir = $(PWD)/../lib
objdir = $(PWD)/../obj
moddir = $(bindir)
srcdir = $(PWD)

# Dependencies
object_spline = \
	SplineAkima.o \
	SplineBessel.o \
	SplineBiCubic.o \
	SplineBilinear.o \
	SplineBiQuintic.o \
	SplineConstant.o \
	SplineCubicBase.o \
	SplineCubic.o \
	SplineHermite.o \
	SplineLinear.o \
	SplinePchip.o \
	SplineQuinticBase.o \
	SplineQuintic.o \
	SplinesBivariate.o \
	Splines.o \
	SplinesCinterface.o \
	SplineSet.o \
	SplineSetGC.o \
	SplinesUnivariate.o \
	SplineVec.o \
	external_spline.o

object_iir = \
	Biquad.o \
	Butterworth.o \
	Cascade.o \
	ChebyshevI.o \
	ChebyshevII.o \
	Custom.o \
	PoleFilter.o \
	RBJ.o \
	external_iir.o \
	module_iirfilt.o

object_filesystem = \
	common.o \
	get_path.o \
	filesystem.o \
	f_filesystem.o \
	f_find.o \
	f_get_path.o \
	f2c.o

object_hdbscan = \
	hdbscan.o \
	cluster.o \
	hdbscanAlgorithm.o \
	hdbscanConstraint.o \
	outlierScore.o \
	undirectedGraph.o \
	hdbscanParameters.o \
	hdbscanResult.o \
	hdbscanRunner.o \
	bitSet.o \
	external_hdbscan.o

object_libflit = \
	module_error.o \
	module_string.o \
	module_constants.o \
	module_array.o \
	module_array_operation.o \
	module_array_extension.o \
	module_sort.o \
	module_calculus.o \
	module_laplacefilt.o \
	module_specialfunc.o \
	module_utility.o \
	module_heap.o \
	module_date_time.o \
	module_unique.o \
	external_random.o \
	module_random.o \
	$(object_filesystem) \
	module_filedir.o \
	module_io.o \
	module_linear_algebra.o \
	module_taper.o \
	$(object_spline) \
	external_mba.o \
	module_interp.o \
	module_fit.o \
	module_random_extra.o \
	external_convexhull.o \
	module_readpar.o \
	module_transform.o \
	module_balancefilt.o \
	module_gaussfilt.o \
	module_statistics.o \
	module_medianfilt.o \
	module_meanfilt.o \
	module_lowessfilt.o \
	module_geometry.o \
	$(object_iir) \
	module_fourierfilt.o \
	module_domain_decomposition.o \
	module_mpicomm.o \
	module_mpicomm_group.o \
	module_andffilt.o \
	module_dipfilt.o \
	module_tvfilt.o \
	module_spectrum.o \
	$(object_hdbscan) \
	module_clustering.o \
	module_libflit.o

object = $(object_libflit)
obj = $(addprefix $(objdir)/, $(object))

# targets
lib = $(bindir)/libflit.a

all: mkdir $(lib) copy clean_genmod

# Compile options
include $(PWD)/Makefile.in

inc = $(base_inc) \
	-I$(PWD) \
	-I$(PWD)/../third_party/pcg_random \
	-I$(PWD)/../third_party/spline \
	-I$(PWD)/../third_party/multilevel_bspline \
	-I$(PWD)/../third_party/iir \
	-I$(PWD)/../third_party/filesystem \
	-I$(PWD)/../third_party/convexhull \
	-I$(PWD)/../third_party/hdbscan/Hdbscan \
	-I$(PWD)/../third_party/hdbscan/HdbscanStar \
	-I$(PWD)/../third_party/hdbscan/Runner \
	-I$(PWD)/../third_party/hdbscan/Utils
lflags = $(base_lflags)
fflags = $(base_fflags) #-O0 -g -check bounds
cflags = $(base_cflags)
cxxflags = $(base_cxxflags)

# Compile
$(objdir)/%.o : $(srcdir)/%.c
	$(cc) -o $@ -c $(cflag) $(inc) $<

$(objdir)/%.o : $(srcdir)/%.cpp
	$(cxx) -o $@ -c $(cxxflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/pcg_random/%.cpp
	$(cxx) -o $@ -c $(cxxflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/spline/%.cc
	$(cxx) -o $@ -c $(cxxflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/iir/%.cpp
	$(cxx) -o $@ -c $(cxxflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/filesystem/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/filesystem/%.c
	$(cc) -o $@ -c $(cflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/filesystem/%.cpp
	$(cxx) -o $@ -c $(cxxflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/hdbscan/*/%.cpp
	$(cxx) -o $@ -c $(cxxflags) $(inc) $<

# link
$(lib): $(obj)
	$(ll) rcs $@ $^

# make directory
mkdir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod $(moddir)/*.smod $(lib)

# copy macro
copy:
	-@cp -p $(srcdir)/libflit.macro $(bindir)

# clean auto-generated f90 files
clean_genmod:
	-@rm -rf $(moddir)/*_genmod.f90
