
# Paths
bindir = $(PWD)/../lib
objdir = $(PWD)/../obj
moddir = $(bindir)
srcdir = $(PWD)

# Dependencies
object_kdtree = \
	kdtree_mod.o

object_mba = \
	mba_mod.o

object_iir = \
	iir_design_mod.o \
	iir_filter_mod.o

object_filesystem = \
	common.o \
	get_path.o \
	filesystem.o \
	f_filesystem.o \
	f_find.o \
	f_get_path.o \
	f2c.o

object_hdbscan = \
	hdbscan_mod.o

object_convexhull = \
	convex_hull_kinds.o \
	convex_hull_2d.o \
	convex_hull_3d.o

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
	$(object_kdtree) \
	$(object_mba) \
	module_interp.o \
	module_fit.o \
	module_random_extra.o \
	module_readpar.o \
	module_transform.o \
	module_balancefilt.o \
	module_gaussfilt.o \
	module_statistics.o \
	module_medianfilt.o \
	module_meanfilt.o \
	$(object_iir) \
	module_iirfilt.o \
	module_lowessfilt.o \
	$(object_convexhull) \
	module_geometry.o \
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
	-I$(PWD)/../third_party/mba \
	-I$(PWD)/../third_party/iir \
	-I$(PWD)/../third_party/kdtree \
	-I$(PWD)/../third_party/filesystem \
	-I$(PWD)/../third_party/convexhull
lflags = $(base_lflags)
fflags = $(base_fflags) #-O0 -g -check bounds -fpe:0 -init:snan
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

$(objdir)/%.o : $(srcdir)/../third_party/kdtree/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/mba/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/iir/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/filesystem/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/filesystem/%.c
	$(cc) -o $@ -c $(cflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/filesystem/%.cpp
	$(cxx) -o $@ -c $(cxxflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/hdbscan/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : $(srcdir)/../third_party/convexhull/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

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
