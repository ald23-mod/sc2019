"""M345SC Homework 2, part 1
Anas Lasri Doukkali, CID: 01209387
"""

#==============================PART 1==========================================#


import numpy as np                                #Importing Numpy

#-----------------------------PART 1.1-----------------------------------------#

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

    Discussion: I will first discuss the implementation and the idea behind the
    algorithm before moving on to analyzing the time complexity of the code itself.

    a) Implementation: The idea will be really simple. I will start in the first
    day by completing the task that have an empty dependency list. These (N-M)
    tasks will then be indexed with the approriate completion date, 0. We now do
    the following, remove the tasks previously completed from any dependency list
    they might be in and complete the tasks that will now have an empty dependency
    list. This process will be repeated until L is emptied of all tasks and all
    all tasks are completed which we keep track of using our task_status list.

    b) Time Complexity: I will now discuss the time complexity of my code and its
    asymptotic behaviour.
    """
    S = np.zeros(len(L))                         # Initializing empty S array
    day_counter = 0                              # Day counter used to update S
    task_status = np.zeros(len(L))               # 0 if not done, 1 otherwise
    while sum(task_status) < len(L):             # While there are undone tasks
        for i in range(len(L)):                  # Looping through dependency list
            if len(L[i]) == 0 and task_status[i] == 0:
                task_status[i] = 1               # Do task if dependency list is empty and still not done
                S[i] = day_counter               # Add the day count to the task done
        for j in range(len(L)):
            if task_status[j] == 1:              # Find tasks previously done in dependency list
                for j2 in range(len(L)):
                    if j in L[j2]:
                        L[j2].remove(j)          # Remove tasks previously done from the dependency list
        day_counter += 1                         # Update the day count
    return S                                     # Return the schedule

#-------------------------------Part 1.2.1-------------------------------------#

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

    L1 = list(np.arange(len(A)))                # Assumes nodes are numbered from 0 to N-1
    L2 = [0 for l in L1]                        # Initializing all nodes as unvisited
    L4 = [[] for l in L1] #paths                # Initializing a path container

    Q=[]                                        # Initializing Queue as empty list
    Q.append(J1)                                # Adding the source node to the queue
    L2[J1]=1                                    # Mark source node as visited
    L4[J1]=[J1]                                 # Path for source
    while len(Q)>0:                             # Stopping algorithm once queue is empty
        x = Q.pop(0)                            # Remove node from front of queue
        print("***x=",x,' ***')
        for i in range(len(A[x])):              # Loop through the adjacency matrix
            v = A[x][i][0]                      # Get the first element of the tuples(node)
            if L2[v]==0 and a0*A[x][i][1] >= amin: # Condition from signal traffic
                Q.append(v)                     # Add unexplored neighbors to back of queue
                L2[v]=1
                #L3[v]=1+L3[x]
                L4[v].extend(L4[x])             # Add path to node x and node v to path
                L4[v].append(v)
            print("v=",v)
            print("Q=",Q)
    L5 = L4[J2]                                 # Printing the correct path

    return L5

#--------------------------------Part 1.2.2------------------------------------#

def a0min(A,amin,J1,J2):
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

    #Initialize dictionaries
    dinit = -1                                # Initial "distance" value
    Edict = {}                                # Explored nodes dict
    Udict = {}                                # Unexplored nodes dict

    for n in range(len(A)):
        Udict[n] = dinit
    Udict[J1] = 0

    #Main Search
    while len(Udict) > 0:
        #Find node with min d in Udict and move to Edict
        dmin  = dinit
        for n,w in Udict.items():
            if w > dmin:
                dmin = w
                nmin = n
        Edict[nmin] = Udict.pop(nmin)
        print("moved node", nmin)

        # Update provisional distances for unexplored neighbours of nmin
        dcomp = 0
        for i in range(len(A[nmin])):
            n = A[nmin][i][0]
            w = A[nmin][i][1]
            if n in Udict:
                if dcomp < w:
                    dcomp = w
                if dcomp < Udict[n] or Udict[n] == -1:
                    Udict[n] = dcomp

    a0_min = amin/Edict[J2]

    path = findPath(A,amin/Edict[J2],amin,J1,J2)

    output = a0_min, path

    if len(path) == 0:
        output = -1,[]


    return output




#if __name__=='__main__':
    #add code here if/as desired
    #L=None #modify as needed
