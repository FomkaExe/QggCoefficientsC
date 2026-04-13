rm -f lsq
rm -rf output
mkdir output
gfortran -std=legacy -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow src/Least_Square3.f -o lsq
