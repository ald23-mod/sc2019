"""M345SC Homework 3, part 1
Anas Lasri Doukkali, CID:01209387
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.signal import hann
import scipy.sparse as sp
import scipy

def nwave(alpha,beta,Nx=256,Nt=801,T=200,display=False):
    """
    Question 1.1
    Simulate nonlinear wave model

    Input:
    alpha, beta: complex model parameters
    Nx: Number of grid points in x
    Nt: Number of time steps
    T: Timespan for simulation is [0,T]
    Display: Function creates contour plot of |g| when true

    Output:
    g: Complex Nt x Nx array containing solution
    """

    #generate grid
    L = 100
    x = np.linspace(0,L,Nx+1)
    x = x[:-1]

    def RHS(f,t,alpha,beta):
        """Computes dg/dt for model eqn.,
        f[:N] = Real(g), f[N:] = Imag(g)
        Called by odeint below
        """
        g = f[:Nx]+1j*f[Nx:]
        #add code here
        #d2g=?
        c = np.fft.fft(g)/Nx
        n = np.arange(-Nx/2,Nx/2)
        n = np.fft.fftshift(n)
        k = 2*np.pi*n/L
        d2g = Nx*np.fft.ifft(-1*k**2*c)
        #-----------
        dgdt = alpha*d2g + g -beta*g*g*g.conj()
        df = np.zeros(2*Nx)
        df[:Nx] = dgdt.real
        df[Nx:] = dgdt.imag
        return df

    #set initial condition
    g0 = np.random.rand(Nx)*0.1*hann(Nx)
    f0=np.zeros(2*Nx)
    f0[:Nx]=g0
    t = np.linspace(0,T,Nt)

    #compute solution
    f = odeint(RHS,f0,t,args=(alpha,beta))
    g = f[:,:Nx] + 1j*f[:,Nx:]

    if display:
        plt.figure()
        plt.contour(x,t,g.real)
        plt.xlabel('x')
        plt.ylabel('t')
        plt.title('Contours of Real(g)')


    return g



def analyze():
    """
    Question 1.2
    Add input/output as needed

    Discussion: Add discussion here
    """

    return None #modify as needed


def wavediff(alpha,a,b,c,beta=0):
    """
    Question 1.3
    Add input/output as needed

    Discussion: Add discussion here
    """

    #----------------------------------------------------
    Nx = 256

    # Calculating g at t = 100 from nwave()
    g = nwave(1-1j,1+2j,Nx,801,100)[800,:]

    # Defining h
    x = np.linspace(0,100,Nx+1)
    x = x[:-1]
    h = x[1] - x[0]

    k = 2*np.pi*4/100
    g = np.sin(k*x)

    hfac = 1/h
    hfac_2 = 1/(2*h)
    hfac_4 = 1/(4*h)
    hfac_6 = 1/(6*h)

    # Defining the row vector b
    b_1 = np.zeros_like(g)
    b_1[1:Nx-1] = hfac_2*a*(g[2:Nx]-g[0:Nx-2])


    b_2 = np.zeros_like(g)
    b_2[2:Nx-2] = hfac_4*b*(g[4:Nx]-g[0:Nx-4])
    b_2[1] = hfac_4*b*(g[3]-g[Nx-1])
    b_2[Nx-2] = hfac_4*b*(g[0]-g[Nx-4])


    b_3 = np.zeros_like(g)
    b_3[3:Nx-3] = hfac_6*c*(g[6:Nx]-g[0:Nx-6])
    b_3[1] = hfac_6*c*(g[4]-g[Nx-2])
    b_3[2] = hfac_6*c*(g[5]-g[Nx-1])
    b_3[Nx-3] = hfac_6*c*(g[0]-g[Nx-6])
    b_3[Nx-2] = hfac_6*c*(g[1]-g[Nx-5])


    db = np.array(b_1) + np.array(b_2) + np.array(b_3)
    #print(db)
    db[0] = hfac*((-17/6)*g[0] + (3/2)*g[1] + (3/2)*g[2] + (-1/6)*g[3])
    db[Nx-1] = -hfac*((-17/6)*g[Nx-1] + (3/2)*g[Nx-2] + (3/2)*g[Nx-3] + (-1/6)*g[Nx-4])


    # Defining the matrix A
    second_row = np.ones(Nx)
    first_row = alpha*np.ones(Nx)
    first_row[0] = 0
    first_row[1] = 3
    third_row = alpha*np.ones(Nx)
    third_row[Nx-1] = 0
    third_row[Nx-2] = 3
    A = np.asarray([first_row,second_row,third_row])


    dg = scipy.linalg.solve_banded((1,1),A,db)

    plt.figure()
    plt.plot(x, k*np.cos(k*x))

    plt.plot(x, dg)
    plt.show()



    return None #modify as needed

if __name__=='__main__':
    x=None
    #Add code here to call functions above and
    #generate figures you are submitting
