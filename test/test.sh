
\rm -rf program_name obj test
\rm -rf test*.txt select*.txt cluster*.txt
\rm -rf exec* *.bin *.png

exit
#
# The following tests only cover some of the functionalities provided by FLIT.
#

# array
touch program_name
mod="array"
echo "program_name="$mod > program_name
make clean
make
./exec_$mod

# filter
touch program_name
mod="filter"
echo "program_name="$mod > program_name
make clean
make
mpirun -np 1 ./exec_$mod
mpirun -np 24 ./exec_$mod

# random
touch program_name
mod="random"
echo "program_name="$mod > program_name
make clean
make
./exec_$mod

# taper
touch program_name
mod="taper"
echo "program_name="$mod > program_name
make clean
make
./exec_$mod

# readpar
touch program_name
mod="readpar"
echo "program_name="$mod > program_name
make clean
make
./exec_$mod \
	par1=0.12345 par2=0:0,10:2.0 par3=1.2,2.3,3.4 \
	cpar1=0.12345+0.1i cpar2=0~3:0.0,10:2.0+3.0i cpar3=1.2,2.3+1.0i,3.4-4.0i

# transform
touch program_name
mod="transform"
echo "program_name="$mod > program_name
make clean
make
./exec_$mod

# interp
touch program_name
mod="interp"
echo "program_name="$mod > program_name
make clean
make
./exec_$mod

# filedir
touch program_name
mod="filedir"
echo "program_name="$mod > program_name
make clean
make
./exec_$mod

# geometry
touch program_name
mod="geometry"
echo "program_name="$mod > program_name
make clean
make
./exec_$mod
