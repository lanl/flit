
# paths
flitdir = $(PWD)/..
bindir = $(PWD)
objdir = ./obj
moddir = ./obj

# dependencies
include $(PWD)/program_name
object = test_$(program_name).o
obj = $(addprefix $(objdir)/, $(object))

# targets
exec = $(PWD)/exec_$(program_name)

all: makedir $(exec)

# options
include $(flitdir)/src/Makefile.in
inc = $(base_inc) \
	-I$(flitdir)/lib
lflags = $(flitdir)/lib/libflit.a \
	$(base_lflags)
fflags = $(base_fflags) #-O0 -g -check bounds
cflags = $(base_cflags)
cxxflags = $(base_cxxflags)

# compile
$(objdir)/%.o : ./%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# link
$(exec) : $(obj)
	$(fc) -o $@ $^ $(lflags) $(inc)

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod
