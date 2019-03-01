"""M345SC Homework 2, part 2
Your name and CID here
"""
import numpy as np                                  # Numpy
import networkx as nx                               # Networkx
from scipy.integrate import odeint                  # Scipy
import matplotlib.pyplot as plt                     # Matplotlib



def model1(G,x=0,params=(50,80,105,71,1,0),tf=6,Nt=400,display=False):
    """
    Question 2.1
    Simulate model with tau=0

    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    tf,Nt: Solutions Nt time steps from t=0 to t=tf (see code below)
    display: A plot of S(t) for the infected node is generated when true

    x: node which is initially infected

    Output:
    S: Array containing S(t) for infected node
    """
    a,theta1,theta2,g,k,tau=params                      # Reading parameters
    tarray = np.linspace(0,tf,Nt+1)                     # Defining the time vector
    S = np.zeros(Nt+1)                                  # Initializing S

    #Function in which the equations will be defined to later be called
    #by an ode solver on Scipy.
    def infection(x,t,params=(50,80,105,71,1,0)):
        a,theta1,theta2,g,k,tau=params
        #Initialize the states
        v = x[0]
        i = x[1]
        s = x[2]

        theta = theta1+theta2*(1-np.sin(2*np.pi*t))     # Precalculate theta
        #We know that the RHS sum will be zero since the flux matrix is zero
        #due to the fact that tau is zero
        dvdt = k*(1 - v) - theta*(s*v)
        didt = theta*s*v - (k + a)*i
        dsdt = a*i - (g + k)*s

        return [dvdt,didt,dsdt]

    x0 = [0.1,0.05,0.05]                                # Initial state conditions
    x = odeint(infection,x0,tarray)                     # Using ode solver

    S = x[:,2]                                          # Retrieving S from entire solution

    if display == True:
        plt.scatter(tarray,S)                           # Scatter plot required
        plt.ylim(-0.005,0.07)
    return S



def modelN(G,x=0,params=(50,80,105,71,1,0.01),tf=6,Nt=400,display=False):
    """
    Question 2.1
    Simulate model with tau=0

    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    tf,Nt: Solutions Nt time steps from t=0 to t=tf (see code below)
    display: A plot of S(t) for the infected node is generated when true

    x: node which is initially infected

    Output:
    Smean,Svar: Array containing mean and variance of S across network nodes at
                each time step.
    """
    a,theta0,theta1,g,k,tau=params              # Reading in parameters
    tarray = np.linspace(0,tf,Nt+1)             # Defining time array
    Smean = np.zeros(Nt+1)                      # Initializing Smean container
    Svar = np.zeros(Nt+1)                       # Initializing Svar container
    N = nx.number_of_nodes(G)                   # Defining number of nodes
    # Initializing y0
    y_S = np.zeros(N)
    y_I = np.zeros(N)
    y_V = np.ones(N)
    y_S[x] = 0.05                               # Spreader ratio initial for infected node
    y_I[x] = 0.05                               # Infected ratio initial for infected node
    y_V[x] = 0.1                                # Vulnerable ratio initial for infected node

    # Concatenating the separate initial condition vectors
    z = np.concatenate([y_S,y_I])
    y0 = np.concatenate([z,y_V])
    #print(y0)


    #-----------------------Defining the adjacency matrix ---------------------#
    A = nx.adjacency_matrix(G).toarray()                   # Adjacency matrix
    q_i = np.array(A.sum(axis=0))                          # Getting the node degrees from Adjacency matrix
    q_im = np.array([q_i,]*N).T                            # Extending previous matrix to square matrix
    sum_kj = np.dot(q_i,A)                                 # Calculating the denominator of the flux matrix
    flux_ij = np.array(tau*np.divide(np.multiply(q_im,A),sum_kj)) # Calculating the flux matrix
    flux_ij[np.isnan(flux_ij)] = 0                         # Nan condition for non connected nodes
    #print(np.sum(flux_ij,axis=0))                         # Can be used to make sure flux matrix is the correct one

    def RHS(y,t):
        """Compute RHS of model at time t
        input: y should be a 3N x 1 array containing with
        S,I,V corresponding to
        S on nodes 0 to N-1, I on nodes 0 to N-1, and
        V on nodes 0 to N-1, respectively.
        output: dy: also a 3N x 1 array corresponding to dy/dt

        Discussion: add discussion here
        """
        #a,theta0,theta1,g,k,tau=params         # Input parameters
        S = np.array(y[:N])                     # Reading initial Spreaders array
        I = np.array(y[N:2*N])                  # Reading initial Infected array
        V = np.array(y[2*N:3*N])                # Reading initial vulnerable array
        dy = np.zeros(3*N)                      # Initializing dy

        #Precalculating
        theta = theta0 + theta1*(1 - np.sin(2*np.pi*t))

        #Setting up the equations
        dy[:N] = a*I - (g + k)*S + np.dot(flux_ij,S) - tau*S                         # S_i
        dy[N:2*N] = theta*np.multiply(S,V) - (k + a)*I + np.dot(flux_ij,I) - tau*I   # I_i
        dy[2*N:3*N] = k - k*V - theta*np.multiply(S,V)  + np.dot(flux_ij,V) - tau*V  # V_i

        return dy

    y = odeint(RHS,y0,tarray)                   # Using odeint scipy ode solver
    S = y[:,:N]                                 # Retrieve S from the solution matrix
    Smean = np.mean(S,axis=1)                   # Calculate the mean of S
    Svar = np.var(S,axis=1)                     # Calculate the variance of Svar


    # If display=True, plot Smean and Svar against time
    if display == True:
        plt.figure()
        plt.scatter(tarray,Smean,s=10,c='blue')
        plt.ylim(0,1.15*max(Smean))
        plt.xlabel('t')
        plt.ylabel('<S(t)>')
        plt.title('Anas Lasri, mean of S(t)')
        plt.savefig('CW2p1.png', dpi = 500)

        plt.figure()
        plt.scatter(tarray,Svar,s=10,c='magenta')
        plt.ylim(0,1.15*max(Svar))
        plt.xlabel('t')
        plt.ylabel('$<(S(t)-<S(t)>)^2>$')
        plt.title('Anas Lasri, variance of S(t)')
        plt.savefig('CW2p2.png', dpi = 500)





    return Smean, Svar

def modelN_modified(G,x=0,params=(0,0,0,0,0,0),tf=6,Nt=400,display=False):
    """
    Question 2.1
    Simulate model with tau=0

    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    tf,Nt: Solutions Nt time steps from t=0 to t=tf (see code below)
    display: A plot of S(t) for the infected node is generated when true

    x: node which is initially infected

    Output:
    Smean,Svar: Array containing mean and variance of S across network nodes at
                each time step.
    """
    a,theta0,theta1,g,k,tau=params
    tarray = np.linspace(0,tf,Nt+1)
    Smean = np.zeros(Nt+1)
    Svar = np.zeros(Nt+1)
    N = nx.number_of_nodes(G)           # Defining number of nodes
    y_S = np.zeros(N)
    y_I = np.zeros(N)
    y_V = np.ones(N)
    y_S[x] = 0.05
    y_I[x] = 0.05
    y_V[x] = 0.1

    z = np.concatenate([y_S,y_I])
    y0 = np.concatenate([z,y_V])


    #-------------------- Defining the adjacency matrix --------------------
    A = nx.adjacency_matrix(G).toarray()
    q_i = np.array(A.sum(axis=0))
    q_im = np.array([q_i,]*N).T
    sum_kj = np.dot(q_i,A)
    flux_ij = np.array(tau*np.divide(np.multiply(q_im,A),sum_kj))
    flux_ij[np.isnan(flux_ij)] = 0



    def RHS(y,t):
        """Compute RHS of model at time t
        input: y should be a 3N x 1 array containing with
        S,I,V corresponding to
        S on nodes 0 to N-1, I on nodes 0 to N-1, and
        V on nodes 0 to N-1, respectively.
        output: dy: also a 3N x 1 array corresponding to dy/dt

        Discussion: add discussion here
        """
        #a,theta0,theta1,g,k,tau=params      # Input parameters
        #N = nx.number_of_nodes(G)           # Defining number of nodes
        S = np.array(y[:N])
        I = np.array(y[N:2*N])
        V = np.array(y[2*N:3*N])
        dy = np.zeros(3*N)                  # Initializing dy

        theta = theta0 + theta1*(1 - np.sin(2*np.pi*t))

        dy[:N] = a*I - (g + k)*S + np.dot(flux_ij,S) - tau*S
        dy[N:2*N] = theta*np.multiply(S,V) - (k + a)*I + np.dot(flux_ij,I) - tau*I
        dy[2*N:3*N] = k - k*V - theta*np.multiply(S,V)  + np.dot(flux_ij,V) - tau*V


        return dy

    y = odeint(RHS,y0,tarray)
    S = y[:,:N]
    I = y[:,N:2*N]
    V = y[:,2*N:3*N]
    Smean = np.mean(S,axis=1)
    Svar = np.var(S,axis=1)
    Imean = np.mean(I,axis=1)
    Ivar = np.var(I,axis=1)
    Vmean = np.mean(V,axis=1)
    Vvar = np.var(V,axis=1)



    return Smean, Svar,Imean, Ivar, Vmean, Vvar


def diffusion(G,tf=100,Nt=400,display=False):
    """Analyze similarities and differences
    between simplified infection model and linear diffusion on
    Barabasi-Albert networks.
    Modify input and output as needed.

    Discussion: add discussion here
    """
    theta_0=  np.linspace(0,30,5)
    tau_0 = np.linspace(0,1,5)
    Smean_m = np.zeros((len(tau_0),Nt+1))
    Svar_m = np.zeros((len(tau_0),Nt+1))
    Vmean_m = np.zeros((len(tau_0),Nt+1))
    Vvar_m = np.zeros((len(tau_0),Nt+1))
    Imean_m = np.zeros((len(tau_0),Nt+1))
    Ivar_m = np.zeros((len(tau_0),Nt+1))

    tarray = np.linspace(0,tf,Nt+1)
    for i in range(len(tau_0)):
        Smean_m[i][:], Svar_m[i][:], Imean_m[i][:], Ivar_m[i][:], Vmean_m[i][:], Vvar_m[i][:] = modelN_modified(G,x=0,params=(0,theta_0[i],0,0,0,tau_0[i]),tf=100,Nt=400,display=False)

    plt.figure()
    for i in range(len(tau_0)):
        plt.hold(True)
        #plt.plot(tarray,Svar_m[i][:])
        plt.plot(tarray,Ivar_m[i][:])
        plt.plot(tarray,Vvar_m[i][:])
    plt.show()

    plt.figure()
    for i in range(len(tau_0)):
        plt.hold(True)
        #plt.plot(tarray,Svar_m[i][:])
        plt.plot(tarray,Imean_m[i][:])
        plt.plot(tarray,Vmean_m[i][:])
    plt.show()






    return None #modify as needed


if __name__=='__main__':
    #add code here to call diffusion and generate figures equivalent
    #to those you are submitting
    G=None #modify as needed
