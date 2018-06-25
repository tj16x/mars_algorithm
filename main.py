# Author:     Thomas Oliver Jelly.
# Location:   University of Glasgow.
# Purpose:    Multi-scale Anisotropic Rough Surface (MARS) algorithm.
# Contact:    thomas.jelly@glasgow.ac.uk

import mars
import numpy as np

# Test imports from hereon
import time
import scipy as sp


"""
Input Variables:
n     : number of streamwise points in the correlation function
m     : number of spanwise   points in the correlation function
N     : number streamwise points on the height map 
M     : number spanwise points on the height map 
dx    : sampling interval in streamwise direction
dy    : sampling interval in spanwise direction
cutoff: cutoff for correlation function
fout  : heightmap filename
"""

n=8
m=n
N=128
M=N
dx=1
dy=1
cutoff=1e-5
fout='heightmap.dat'

# Create an instance of the surface class
s= mars.surface(n,m,M,N,cutoff,dx,dy)

# Step 1: Specify ACF
acf= s.acf()

# Step 2: Assemble & solve nonlinear system of equations
guess= s.f0()

time1 = time.clock()
alpha= s.ncgm(guess)
time2 = time.clock()


# Step 3: Generate a random number matrix
rand= s.eta()

# Step 4: Generate the heightmap
hmap= s.heightmap(alpha,rand)

# Step 5: Save the surface
s.save(fout)
