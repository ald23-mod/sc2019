"""M345SC Homework 2, part 1
Anas Lasri Doukkali, CID: 01209387
"""
import numpy as np
from collections import defaultdict

def scheduler(L):
    """
    Question 1.1
    Schedule tasks using dependency list provided as input

    Input:
    L: Dependency list for tasks. L contains N sub-lists, and L[i] is a sub-list
    containing integers (the sub-list my also be empty). An integer, j, in this
    sub-list indicates that task j must be completed before task i can be started.

    Output:
    S: A list of integers corresponding to the schedule of tasks. S[i] indicates
    the day on which task i should be carried out. Days are numbered starting
    from 0.

    Discussion: Add analysis here
    """
    S= [0]*len(L)
    day_counter = 0
    task_status = [0] * len(L)
    while sum(task_status) < len(L):
        for i in range(len(L)):
            if len(L[i]) == 0 and task_status[i] == 0:
                task_status[i] = 1
                S[i] = day_counter
        for j in range(len(L)):
            if task_status[j] == 1:
                for j2 in range(len(L)):
                    if j in L[j2]:
                        L[j2].remove(j)
        day_counter += 1
    return S


def findPath(A,a0,amin,J1,J2):
    """
    Question 1.2 i)
    Search for feasible path for successful propagation of signal
    from node J1 to J2

    Input:
    A: Adjacency list for graph. A[i] is a sub-list containing two-element tuples (the
    sub-list my also be empty) of the form (j,Lij). The integer, j, indicates that there is a link
    between nodes i and j and Lij is the loss parameter for the link.

    a0: Initial amplitude of signal at node J1

    amin: If a>=amin when the signal reaches a junction, it is boosted to a0.
    Otherwise, the signal is discarded and has not successfully
    reached the junction.

    J1: Signal starts at node J1 with amplitude, a0
    J2: Function should determine if the signal can successfully reach node J2 from node J1

    Output:
    L: A list of integers corresponding to a feasible path from J1 to J2.

    Discussion: Add analysis here
    """

    L1 = list(np.arange(len(A)))      #Assumes nodes are numbered from 0 to N-1
    L2 = [0 for l in L1]
    L3 = [-1000 for l in L1]
    L4 = [[] for l in L1] #paths


    Q=[]
    Q.append(J1)
    L2[J1]=1
    L3[J1]=0
    L4[J1]=[J1]
    while len(Q)>0:
        x = Q.pop(0) #remove node from front of queue
        print("***x=",x,' ***')
        for i in range(len(A[x])):
            v = A[x][i][0]
            if L2[v]==0 and a0*A[x][i][1] >= amin:
                Q.append(v) #add unexplored neighbors to back of queue
                L2[v]=1
                L3[v]=1+L3[x]
                L4[v].extend(L4[x]) #Add path to node x and node v to path
                L4[v].append(v)     #for node v


            print("v=",v)
            print("Q=",Q)

    L5 = L4[J2]

    return L5



def a0min(A,amin,J1,J2,a0):
    """
    Question 1.2 ii)
    Find minimum initial amplitude needed for signal to be able to
    successfully propagate from node J1 to J2 in network (defined by adjacency list, A)

    Input:
    A: Adjacency list for graph. A[i] is a sub-list containing two-element tuples (the
    sub-list my also be empty) of the form (j,Lij). The integer, j, indicates that there is a link
    between nodes i and j and Lij is the loss parameter for the link.

    amin: Threshold for signal boost
    If a>=amin when the signal reaches a junction, it is boosted to a0.
    Otherwise, the signal is discarded and has not successfully
    reached the junction.

    J1: Signal starts at node J1 with amplitude, a0
    J2: Function should determine min(a0) needed so the signal can successfully
    reach node J2 from node J1

    Output:
    (a0min,L) a two element tuple containing:
    a0min: minimum initial amplitude needed for signal to successfully reach J2 from J1
    L: A list of integers corresponding to a feasible path from J1 to J2 with
    a0=a0min
    If no feasible path exists for any a0, return output as shown below.

    Discussion: Add analysis here
    """
    L1 = list(np.arange(len(A)))      #Assumes nodes are numbered from 0 to N-1
    L2 = [0 for l in L1]
    L3 = [-1000 for l in L1]
    L4 = [[] for l in L1] #paths


    Q=[]
    Q.append(J1)
    L2[J1]=1
    L3[J1]=0
    L4[J1]=[J1]
    while len(Q)>0:
        x = Q.pop(0) #remove node from front of queue
        print("***x=",x,' ***')
        for i in range(len(A[x])):
            v = A[x][i][0]
            if L2[v]==0 and a0*A[x][i][1] >= amin:
                Q.append(v) #add unexplored neighbors to back of queue
                L2[v]=1
                L3[v]=1+L3[x]
                L4[v].extend(L4[x]) #Add path to node x and node v to path
                L4[v].append(v)     #for node v


            print("v=",v)
            print("Q=",Q)

    L5 = L4[J2]


    #output = -1,[] #Modify as needed

    return #output


#if __name__=='__main__':
    #add code here if/as desired
    #L=None #modify as needed
