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
    x0 = [0.1,0.05,0.05]
    x = odeint(infection,x0,tarray)

    S = x[:,2]

    plt.scatter(tarray,S)
    plt.ylim(-0.005,0.07)
    return S

def infection(x,t,params=(50,80,105,71,1,0)):
    a,theta1,theta2,g,k,tau=params
    v = x[0]
    i = x[1]
    s = x[2]

    theta = theta1+theta2*(1-np.sin(2*np.pi*t))

    didt = theta*s*v - (k + a)*i
    dvdt = k*(1 - v) - theta*(s*v)
    dsdt = a*i - (g + k)*s

    return [dvdt,didt,dsdt]


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
    a,theta1,theta2,k,g,tau=params
    tarray = np.linspace(0,tf,Nt+1)
    Smean = np.zeros(Nt+1)
    Svar = np.zeros(Nt+1)

    #Add code here

    def RHS(y,t):
        """Compute RHS of model at time t
        input: y should be a 3N x 1 array containing with
        y[:N],y[N:2*N],y[2*N:3*N] corresponding to
        S on nodes 0 to N-1, I on nodes 0 to N-1, and
        V on nodes 0 to N-1, respectively.
        output: dy: also a 3N x 1 array corresponding to dy/dt

        Discussion: add discussion here
        """

        dy = 0 #modify
        return dy



    return Smean,Svar


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
