import numpy as np
import matplotlib.pyplot as plt

fname = "../heightmap.dat"
with open(fname, 'r') as f:
   header = [next(f) for x in xrange(4)]

Nx = int(header[0].replace("\n","").strip())
Ny = int(header[1].replace("\n","").strip())
Lx = float(header[2].replace("\n","").strip())
Ly = float(header[3].replace("\n","").strip())


data = np.genfromtxt(fname, skip_header=4)

X = data[:,0].reshape(Nx, Ny)
Y = data[:,1].reshape(Nx, Ny)
heightmap = data[:,2].reshape(Nx, Ny)

plt.pcolormesh(X, Y, heightmap, cmap=plt.cm.RdYlBu_r)
plt.axis("tight")
plt.show()
