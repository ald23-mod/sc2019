"""M345SC Homework 3, part 2
Anas Lasri Doukkali, CID:01209387
"""
import numpy as np
import networkx as nx
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def growth1(G,params=(0.02,6,1,0.1,0.1),T=6):
    """
    Question 2.1
    Find maximum possible growth, G=e(t=T)/e(t=0) and corresponding initial
    condition.

    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    T: time for energy growth

    Output:
    G: Maximum growth
    y: 3N-element array containing computed initial condition where N is the
    number of nodes in G and y[:N], y[N:2*N], y[2*N:] correspond to S,I,V

    Discussion: Add discussion here
    """
    a,theta,g,k,tau=params
    N = G.number_of_nodes()

    #Construct flux matrix (use/modify as needed)
    Q = [t[1] for t in G.degree()]

    Pden = np.zeros(N)
    Pden_total = np.zeros(N)
    for j in G.nodes():
        for m in G.adj[j].keys():
            Pden[j] += Q[m]
    Pden = 1/Pden
    Q = tau*np.array(Q)
    F = nx.adjacency_matrix(G).toarray()
    F = F*np.outer(Q,Pden)
    #-------------------------------------
    G=0
    y = np.zeros(3*N)

    #Add code here
    M_1 = scipy.linalg.block_diag(F,F,F)


    one_1 = np.array((-g-k-tau)*np.ones(N))
    two_2 = np.array((-k-a-tau)*np.ones(N))
    three_3 = np.array((k-tau)*np.ones(N))
    alpha_d = np.array(np.diag(np.array(a*np.ones(N))))
    theta_d = np.array(np.diag(np.array(theta*np.ones(N))))
    theta_2d = np.array(np.diag(np.array(-theta*np.ones(N))))

    two_d = np.diag(two_2)
    one_d = np.diag(one_1)
    three_d = np.diag(three_3)
    zeros = np.zeros((N,N))

    M_2 = np.block([[one_d,alpha_d,zeros],[theta_d,two_d,zeros],[theta_2d,zeros,three_d]])
    M = np.array(M_1 + M_2)

    #-----------------------------------------------
    A = scipy.linalg.expm(M*T)
    u,y,vh = np.linalg.svd(A)

    x_t = np.dot(A,vh[0])

    growth_t = np.linalg.norm(x_t)**2
    growth_0 = np.linalg.norm(vh[0])

    max_growth = growth_t/growth_0

    return max_growth,vh[0]

def growth2(G,params=(0.02,6,1,0.1,0.1),T=6):
    """
    Question 2.2
    Find maximum possible growth, G=sum(Ii^2)/e(t=0) and corresponding initial
    condition.

    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    T: time for energy growth

    Output:
    G: Maximum growth
    y: 3N-element array containing computed initial condition where N is the
    number of nodes in G and y[:N], y[N:2*N], y[2*N:] correspond to S,I,V

    Discussion: Add discussion here
    """
    a,theta,g,k,tau=params
    N = G.number_of_nodes()

    #Construct flux matrix (use/modify as needed)
    Q = [t[1] for t in G.degree()]

    Pden = np.zeros(N)
    Pden_total = np.zeros(N)
    for j in G.nodes():
        for m in G.adj[j].keys():
            Pden[j] += Q[m]
    Pden = 1/Pden
    Q = tau*np.array(Q)
    F = nx.adjacency_matrix(G).toarray()
    F = F*np.outer(Q,Pden)
    #-------------------------------------
    G=0
    y = np.zeros(3*N)

    #Add code here
    M_1 = scipy.linalg.block_diag(F,F,F)


    one_1 = np.array((-g-k-tau)*np.ones(N))
    two_2 = np.array((-k-a-tau)*np.ones(N))
    three_3 = np.array((k-tau)*np.ones(N))
    alpha_d = np.array(np.diag(np.array(a*np.ones(N))))
    theta_d = np.array(np.diag(np.array(theta*np.ones(N))))
    theta_2d = np.array(np.diag(np.array(-theta*np.ones(N))))

    two_d = np.diag(two_2)
    one_d = np.diag(one_1)
    three_d = np.diag(three_3)
    zeros = np.zeros((N,N))

    M_2 = np.block([[one_d,alpha_d,zeros],[theta_d,two_d,zeros],[theta_2d,zeros,three_d]])
    M = np.array(M_1 + M_2)


    #-----------------------------------------------
    A = scipy.linalg.expm(M*T)
    print(A.shape)
    u,y_b,y = np.linalg.svd(A)

    x_t = np.dot(A[N:2*N,:],y[0])

    growth_t = np.linalg.norm(x_t)**2
    growth_0 = np.linalg.norm(y[0])**2

    G = growth_t/growth_0

    return G,y


def growth3(G,params=(2,2.8,1,1.0,0.5),T=6):
    """
    Question 2.3
    Find maximum possible growth, G=sum(Si Vi)/e(t=0)
    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    T: time for energy growth

    Output:
    G: Maximum growth

    Discussion: Add discussion here
    """
    a,theta,g,k,tau=params
    N = G.number_of_nodes()

    #Construct flux matrix (use/modify as needed)
    Q = [t[1] for t in G.degree()]

    Pden = np.zeros(N)
    Pden_total = np.zeros(N)
    for j in G.nodes():
        for m in G.adj[j].keys():
            Pden[j] += Q[m]
    Pden = 1/Pden
    Q = tau*np.array(Q)
    F = nx.adjacency_matrix(G).toarray()
    F = F*np.outer(Q,Pden)
    #-------------------------------------
    G=0
    y = np.zeros(3*N)

    #Add code here
    M_1 = scipy.linalg.block_diag(F,F,F)


    one_1 = np.array((-g-k-tau)*np.ones(N))
    two_2 = np.array((-k-a-tau)*np.ones(N))
    three_3 = np.array((k-tau)*np.ones(N))
    alpha_d = np.array(np.diag(np.array(a*np.ones(N))))
    theta_d = np.array(np.diag(np.array(theta*np.ones(N))))
    theta_2d = np.array(np.diag(np.array(-theta*np.ones(N))))

    two_d = np.diag(two_2)
    one_d = np.diag(one_1)
    three_d = np.diag(three_3)
    zeros = np.zeros((N,N))

    M_2 = np.block([[one_d,alpha_d,zeros],[theta_d,two_d,zeros],[theta_2d,zeros,three_d]])
    M = np.array(M_1 + M_2)


    #-----------------------------------------------
    A_s = scipy.linalg.expm(M*T)[:N,:]
    A_v = scipy.linalg.expm(M*T)[2*N:3*N,:]

    A = np.dot(A_s.T,A_v)

    w,v = np.linalg.eigh(0.5*(A+A.T))

    growth_t = max(w)

    G = growth_t


    return G


def Inew(D):
    """
    Question 2.4

    Input:
    D: N x M array, each column contains I for an N-node network

    Output:
    I: N-element array, approximation to D containing "large-variance"
    behavior

    Discussion: Add discussion here
    """
    D = D
    N,M = D.shape
    I = np.zeros(N)
    D2 = D - np.outer(np.ones((N,1)),D.mean(axis=0))
    U,S,VT = np.linalg.svd(D2.T)
    G = np.dot(U.T,D2.T)
    I = G[0]
    #Add code here
    #fig = plt.figure()
    #ax = Axes3D(fig, elev=70, azim=135)
    #for i in [0,50,100]:
    #    ax.scatter(G[0,i:i+50],G[1,i:i+50],G[2,i:i+50])
    #plt.show()
    return I


if __name__=='__main__':
    G=None
    #add/modify code here if/as desired
    #N,M = 100,5
    #G = nx.barabasi_albert_graph(N,M,seed=1)
