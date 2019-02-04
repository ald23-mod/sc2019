""" Simple bfs code from lecture 7
"""
import networkx as nx

def bfs(G,x):
    """
    Input:
    G: networkx graph
    x: source node

    Output:
    L2: Labels for all nodes in graph, 0=unreachable from source, 1=reachable
    L3: Shortest distance from source to nodes in graph
    """

    L1 = list(G.nodes) #Assumes nodes are numbered from 0 to N-1
    L2 = [0 for l in L1]
    L2 = [0 for l in L1]
    L3 = [-1000 for l in L1]

    Q=[]
    s = 1
    Q.append(s)
    L2[s]=1
    L3[s]=0
    while len(Q)>0:
        x = Q.pop(0)
        for v in G.adj[x].keys():
            if L2[v]==0:
                Q.append(v)
                L2[v]=1
                L3[v]=1+L3[x]
            print("v=",v)
            print("Q=",Q)
    return L2,L3
