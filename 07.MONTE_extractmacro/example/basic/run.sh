#!/bin/bash

# Reset example by removing files 
rm -r material_xs

# Run MONTE_extractmacro.x
../../MONTE_extractmacro.x

# run python script to get mgxs.h5 file. Make sure you have openmc installed
# python3 ../xsfiles2openmc.py
