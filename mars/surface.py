"""
  __  __          _____   _____            _      _____  ____  _____  _____ _______ _    _ __  __  
 |  \/  |   /\   |  __ \ / ____|     /\   | |    / ____|/ __ \|  __ \|_   _|__   __| |  | |  \/  | 
 | \  / |  /  \  | |__) | (___      /  \  | |   | |  __| |  | | |__) | | |    | |  | |__| | \  / | 
 | |\/| | / /\ \ |  _  / \___ \    / /\ \ | |   | | |_ | |  | |  _  /  | |    | |  |  __  | |\/| | 
 | |  | |/ ____ \| | \ \ ____) |  / ____ \| |___| |__| | |__| | | \ \ _| |_   | |  | |  | | |  | | 
 |_|  |_/_/    \_\_|  \_\_____/  /_/    \_\______\_____|\____/|_|  \_\_____|  |_|  |_|  |_|_|  |_| 
                                                                                                   
 Contact: thomas.jelly@glasgow.ac.uk                                                T.O.J., 2018. 
 Purpose: Numerical generation of realistic rough surfaces with specified statistical parameters.

"""

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, root

# Surface class
class surface():
    " Surface class for MARS algorithm "

    # Constructor
    def __init__(self, n=0, m=0, N=0, M=0, c=0, dx=0, dy=0, phi=0):
        """
        Input arguments:
            m (int)
            n (int)
            M (int)
            N (int)
            c (float)
        """
        print "\nWelcome to MARS."
        self.n = n
        self.m = m
        self.N = N
        self.M = M 
        self.c = c
        self.dx = dx
        self.dy = dy
        self.phi = phi
        self.iterations = 0

    # Assemble autocorrelation coefficient function (ACF)
    def acf(self):
        " Assemble [n,m] autocorrelation coefficient function "

        lx = float(self.n-1)
        ly = float(self.m-1)

        exp = lambda x,y,X,Y: np.exp(-2.3*np.sqrt((x/X)**2.0+(y/Y)**2.0 ))
        dx,dy = np.meshgrid(np.arange(0,self.M),np.arange(0,self.N))
        acf = np.zeros([self.N,self.M])
        
        # Evaluate the arbitrary direction of the function 
        #  ->Equation 36 from Bakalos (2003)
        for p in range(self.N):
            for q in range(self.M):
                
                ca = np.cos(self.phi * np.pi/180.0)
                sa = np.sin(self.phi * np.pi/180.0)
                xp =  p*ca + q*sa
                yp = -p*sa + q*ca
                
                acf[p,q] = exp(xp,yp,lx,ly)


        idx = lambda arr,vec: (np.abs(arr-vec)).argmin()+1
        self.n = idx(acf[:,0],self.c)
        self.m = idx(acf[0,:],self.c)

        self.rhs = acf[0:self.n,0:self.m]

        return self.rhs

    # Assemble initial guess for solution to non-linear system of equations
    def f0(self):
        " Assemble [n,m] initial guess for iterative solution "

        d1,d2 = np.meshgrid(np.ones(self.m)*self.m - np.arange(self.m),
                            np.ones(self.n)*self.n - np.arange(self.n) )

        c = np.zeros([self.n,self.m])
        c[:,:] = self.rhs[:,:]/(d1*d2)
        s2 = 1.0/np.sum(c**2.0)

        self.guess = np.sqrt(s2)*c

        return self.guess

    # Assemble the non-linear system of equations
    def f(self, alpha, acf):
        " Assemble [n,m] system of non-linear simultaneous equations "

        alpha = np.reshape(alpha,[self.n,self.m])

        self.fx = np.zeros([self.n,self.m])
        for p in range(self.n):
            for q in range(self.m):
                self.fx[p][q] += np.sum(alpha[0:self.n-p, 0:self.m-q]*alpha[p:self.n, q:self.m],(0,1))

        self.fx -= self.rhs

        return self.fx.flatten()

    # Fletcher-Reeves-Secant variant of the non-linear conjugate gradient algorithm
    def ncgm(self,x0):  
        " Determine coefficients using non-linear conjugate gradient algorithm "

        x= x0.flatten()

        i = 0 
        k = 0
        n = 4
        s0 = 1e-3
        imax = 1024
        jmax = 2
        cg_eps = 1e-4
        sc_eps = 1e-9

        fp = self.f(x0,self.rhs)
        r = -1.0*fp
        d = r 
        delta_new = np.dot(r,d)
        delta_old = delta_new

        print "\nNonlinear conjugate gradient algorithm initialised..."

        tol_outer = (cg_eps**2.0)*delta_old
        while (i < imax) and (delta_new > tol_outer):

            j = 0
            delta_d = np.dot(d,d)
            alpha = -1.0*s0
            fp1 = self.f(x+s0*d,self.rhs)
            eta_prev = np.dot(fp1,d) 
        
            while True:

                eta = np.dot(fp,d)
                alpha *= eta/(eta_prev-eta)
                x += alpha*d
                eta_prev = eta
                j+= 1 

                tol_inner = alpha*alpha*delta_d
                if (j < jmax) or (tol_inner > sc_eps**2.0):
                    break

            fp = self.f(x,self.rhs)
            r = -1.0*fp

            delta_old = delta_new
            delta_new = np.dot(r,r)
            beta = delta_new/delta_old
            k += 1

            if (k==n) or (np.dot(r,d)<0):
                d = r
                k = 0
            i += 1 

            if (i==1 or i%50==0):
                print "\n  iter","NCGM residual"

            if (i==1 or i%10==0):
                print "  %03d" % i, " %.7e" % delta_new

            if (delta_new > 1e+9):
                print "\nWARNING: Nonlinear conjugate gradient algorithm has NOT converged successfully after",i,"iterations."
                print   "WARNING: Coefficients will be determined using inbuilt function fsolve... This may take some time!"
                x = fsolve(self.f,self.f0().flatten(),(self.acf),fprime=self.fjacobian,xtol=1e-7)
                break

        if ( i<imax ) and ( delta_new < tol_outer ):
            print "\nNonlinear conjugate gradient algorithm has converged successfully after",i,"iterations.\n"

        alpha= np.reshape(x,[self.n,self.m])

        self.residual(alpha)

        return alpha
    
    # Krylov approximation for inverse Jacobian
    def krylov(self, tolerance=1e-7):
        print "\nKrylov method initialised...\n"
        
        optionsList = {'xatol':1e-7} #, 'maxiter':100000

        x = root(self.f, self.f0().flatten(), self.acf, method='krylov',  tol=tolerance, callback=self.plot_residual, options=optionsList)
        if x['success']:
            print(x['message'][:-1]+" after " +str(x['nit']) + " iterations.\n")
        else:
            print("WARNING: The Krylov algorithm did NOT converge succesfully.\n")
            print(x['message'])
            return
        alpha = np.reshape(x['x'], [self.n,self.m])
    
        self.residual(alpha)

        return alpha
    
    def plot_residual(self, solution, residual):
        plt.scatter(self.iterations, residual[0], marker='D', edgecolors='k', c='#FF33FF')
        
        plt.xlabel("Iteration")
        plt.ylabel("Residual")
        plt.title("Last Residual: "+ str(residual[0]))
        
        plt.pause(0.05)
        self.iterations += 1
    

    # Assemble the Jacobian for solution to non-linear system of equations
    def fjacobian(self,alpha,acf):
        " Assemble [n*m,n*m] Jacobian "

        alpha = np.reshape(alpha,[self.n,self.m])

        j1 = np.zeros([self.n*self.m,self.n*self.m])
        for p in range(self.n):
            for q in range(self.m):
                for i in range(1,self.n-p+1):
                    for j in range(1,self.m-q+1):
                        r = (p*self.m)+q+1
                        s = (i-1)*self.m+j
                        j1[r-1,s-1] = alpha[i+p-1,j+q-1]

        j2 = np.zeros([self.n*self.m,self.n*self.m])
        for p in range(self.n):
            for q in range(self.m):
                for i in range(1+p,self.n+1):
                    for j in range(1+q,self.m+1):
                        r = (p*self.m)+q+1
                        s = (i-1)*self.m+j
                        j2[r-1,s-1] = alpha[i-p-1,j-q-1]

        self.jacobian = j1+j2

        return self.jacobian

    # Quantify the residual |Ax-b_{sol}|
    def residual(self,alpha):
        " Quantify residual of numerical solution "

        self.res = np.zeros([self.n,self.m])
        for p in range(self.n):
            for q in range(self.m):
                self.res[p][q] += np.sum(alpha[0:self.n-p, 0:self.m-q]*alpha[p:self.n, q:self.m],(0,1))

        self.res -= self.rhs

        l2 = np.sqrt(np.sum(np.absolute(self.res)**2.0))
        print "The l2 norm of the residual is",l2

        return self.res

    # Generate a Gaussian random number matrix
    def eta(self):

        mu = 0.0
        sigma = 1.0
        self.rand = np.random.normal(mu,sigma,[self.N+self.n,self.M+self.m])
    
        return self.rand

    # Generate the heightmap
    def heightmap(self,alpha,rand):
        " Generate heightmap via linear transformation "

        ko = np.mod(np.arange(self.N+self.n),self.N)
        lo = np.mod(np.arange(self.M+self.m),self.M)

        self.hmap = np.zeros([self.N,self.M])
        for i in range(self.N):
            for j in range(self.M):
                for k in range(self.n):
                    for l in range(self.m):
                        self.hmap[i,j] += alpha[k,l]*rand[ko[i+k],lo[j+l]]

        return self.hmap

    # Save the heightmap
    def save(self,fname):
        " Write heightmap to disk "
  
        x,y = np.meshgrid(np.arange(self.N),np.arange(self.M))

        h1 = np.matrix([self.N,self.M]).T
        h2 = np.matrix([self.dx,self.dy]).T
        mat = np.matrix([x.flatten(), y.flatten(), self.hmap.flatten()]).T

        with open(fname,'w') as f:
            np.savetxt(f, h1, fmt=['%-7d'])
            np.savetxt(f, h2, fmt=['%-7.12f'])
            for row in mat:
                np.savetxt(f, row, fmt=['%-7d','%-7d','% 7.12f'])

        print "\nA heightmap called " + fname + " has been saved succesffully.\n"

        return
