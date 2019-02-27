"""M345SC Homework 2, part 2
Your name and CID here
"""
import numpy as np
import networkx as nx
from scipy.integrate import odeint
import matplotlib.pyplot as plt



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
    a,theta1,theta2,g,k,tau=params
    tarray = np.linspace(0,tf,Nt+1)
    S = np.zeros(Nt+1)

    #Add code here
    def infection(x,t,params=(50,80,105,71,1,0)):
        a,theta1,theta2,g,k,tau=params
        v = x[0]
        i = x[1]
        s = x[2]

        theta = theta1+theta2*(1-np.sin(2*np.pi*t))

        dvdt = k*(1 - v) - theta*(s*v)
        didt = theta*s*v - (k + a)*i
        dsdt = a*i - (g + k)*s

        return [dvdt,didt,dsdt]
    x0 = [0.1,0.05,0.05]
    x = odeint(infection,x0,tarray)

    S = x[:,2]

    plt.scatter(tarray,S)
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
    #print(y0)


    #-------------------- Defining the adjacency matrix --------------------
    A = nx.adjacency_matrix(G).toarray()
    print(A)
    q_i = np.array(A.sum(axis=0))
    q_im = np.array([q_i,]*N).T
    sum_kj = np.dot(q_i,A)
    #sum_kjm = np.array([sum_kj,]*N).T
    flux_ij = np.array(tau*np.divide(np.multiply(q_im,A),sum_kj))
    flux_ij[np.isnan(flux_ij)] = 0
    print(flux_ij)
    print(sum_kj)
    #print(tau*np.multiply(q_im,A))
    #print(sum_kjm)
    #print(flux_ij)



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




        dy[:N] = a*I - (g + k)*S + np.dot(flux_ij,S) - np.multiply(np.array(flux_ij.sum(axis=0)),S)
        dy[N:2*N] = theta*np.multiply(S,V) - (k + a)*I + np.dot(flux_ij,I) - np.multiply(np.array(flux_ij.sum(axis=0)),I)
        dy[2*N:3*N] = k - k*V - theta*np.multiply(S,V)  + np.dot(flux_ij,V) - np.multiply(np.array(flux_ij.sum(axis=0)),V)


        return dy

    y = odeint(RHS,y0,tarray)
    S = y[:,:N]
    Smean = np.mean(S,axis=1)
    Svar = np.var(S,axis=1)

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





    return #Smean, Svar


def diffusion(input=(None)):
    """Analyze similarities and differences
    between simplified infection model and linear diffusion on
    Barabasi-Albert networks.
    Modify input and output as needed.

    Discussion: add discussion here
    """


    return None #modify as needed


if __name__=='__main__':
    #add code here to call diffusion and generate figures equivalent
    #to those you are submitting
    G=None #modify as needed
