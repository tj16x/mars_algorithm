# Author:     Thomas Oliver Jelly.
# Location:   University of Glasgow.
# Purpose:    Multi-scale Anisotropic Rough Surface (MARS) algorithm.
# Contact:    thomas.jelly@glasgow.ac.uk

import mars
import numpy as np
import sys

"""
Input Variables:
n     : number of streamwise points in the correlation function
m     : number of spanwise   points in the correlation function
N     : number streamwise points on the height map 
M     : number spanwise points on the height map 
dx    : sampling interval in streamwise direction
dy    : sampling interval in spanwise direction
cutoff: cutoff for correlation function
phi   : rotation angle for anisotropic surfaces
skew  : target skewness for the surface
kurt  : target kurtosis for the surface
fout  : heightmap filename
"""

n=5
m=n
N=128
M=N
dx=1
dy=1
cutoff=1e-5
phi=0
skew=0.3
kurt=3.5
fout='heightmap.dat'


# Create an instance of the surface class
s= mars.surface(n,m,N,M,cutoff,dx,dy,phi)

# Step 1: Specify ACF
acf= s.acf()
s.plot_contour(acf)

# Step 2: Assemble & solve nonlinear system of equations
guess= s.f0()

alpha= s.krylov(method="lgmres")
if len(alpha) == 1 and alpha[0] == -1:
    sys.exit("Error: The program couldn't finish execution.")
# Step 3: Generate a random number matrix
rescaled_skew, rescaled_kurt = s.rescale(alpha, skew, kurt)
if (np.isnan(rescaled_skew)) and (np.isnan(rescaled_kurt)):
    sys.exit("Error: The program can't generate the surface because of the input Skew and Kurtosis.")
rand= s.johnson_eta(rescaled_skew, rescaled_kurt)

# Step 4: Generate the heightmap
hmap= s.heightmap(alpha,rand)

# Step 5: Save the surface
s.save(fout)
