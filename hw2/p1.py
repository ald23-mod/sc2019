"""M345SC Homework 2, part 1
Anas Lasri Doukkali, CID: 01209387
"""
#==============================================================================#
#==============================PART 1==========================================#
#==============================================================================#

import numpy as np                                #Importing Numpy

#------------------------------------------------------------------------------#
#-----------------------------PART 1.1-----------------------------------------#
#------------------------------------------------------------------------------#

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
    asymptotic behaviour. We can see we straight away that our code depends heavily
    on "len(L)", hence, all the loops we have in our code run through this length
    making it thus our N. We start by initializing the day count which is constant
    time, the initialization we will need to take into account is the one for the
    zeros matrices which is done for both S and task_status. This is tricky to estimate
    as it will depend on the computer architecture. For computers with lower Cache(L1, L2...)
    storage there will be abrupt jumps as the matrices will need to be stored in the
    bigger sized memories which are normaly further away from the arithmetic units.
    This will however not affect our leading order or time complexity analysis as we
    will be dealing with higher orders. We start with a while loop, the order for this
    while loop depends on heavily on the structure of the input and will be mainly
    carachterized by the longest dependency chain, hence, if we get a list in which the
    number of tasks with an empty dependency list far exceeds the number of tasks that
    have a non-empty dependency list we then consider that our largest dependency list,
    call it D will be D<<N, however the worst case will be when the longest dependency list
    is approximately N, taking into account that the maximum it could be is N-1.
    Now, inside this while loop we have two separate for loops.
    In the first for loop, which has order N we simply include an if statement which then
    only execute simple tasks if they happen to be true, hence this will not affect the
    leading order or time complexity. The second for loop however has a nested for loop
    which comes from an if statement, hence this for loop will not be excecuted always
    however as we will be considering the worst case scenario we can say that with this
    final for loop our code can have a leading order of O(N^3), however this will not be the
    case regularly as we have if statement that amortisize our time complexity.
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

#------------------------------------------------------------------------------#
#-------------------------------Part 1.2.1-------------------------------------#
#------------------------------------------------------------------------------#

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

    Discussion: I will first discuss the implementation and the idea behind the
    algorithm before moving on to analyzing the time complexity of the code itself.

    a) Implementation: The approach taken for this part of the coursework was the use of a modified
    bfs or dfs to be able to find any feasable path and then print said path. The
    code below will be for bfs however this can easily be modified for dfs by modifying
    the pop argument to be empty. The way this works is as described in lectures.
    We start as always by initialization of some needed lists. L1 will be a list of the nodes,
    while L2 will be a list used to denote whether the nodes have been visited or no.
    Hence if the node i is unvisited then L2[i]=0 and otherwise 1. L4 on the other hand will
    contain the paths which we will be saving as we go along.
    As I mentioned earlier, this code is similar to the one provided in the folder "codes"
    and was only changed for two reasons. The one in lectures uses networkx and we now
    have to write a code that does not. The second reason is that the code from lectures
    only gives us the shortest distance to the destination node. While what we need is the
    path itself. This was easily done by recalling Lab4 where we actually implemented a
    solution to such problem. I choose breadth first search as it is the one that will
    find the shortest path which dfs does not. This in our case was not particularly helpful
    as we were only asked to find a "feasible" path, this is however worth mentioning. ONe
    more addition to the bfs algorithm provided was the condition check in which we now check
    whether the node has enough signal left to be able to make it to the next junction on top
    of whether it was marked as visited or not previously. This is quite important as it eliminates
    all the path in which the signal reaches below minimal levels and hence does not make it across.

    b) Time Complexity: the code only contains two loops. An initial while loop with
    a nested for loop inside. The while loop runs until our queue is empty. We will now
    again need to consider the worst case scenario for this which is when we will have to visit
    all the nodes that we have. Let us say this is order N. All the operation until the
    while loop will not be relevant when taking N large, and hence we will not take them into
    account. Inside this while loop now we have as we said earlier a for loop in which in the
    best case we only execute one operation or in the worst case run four. This loop however
    will do the number of operations mentioned above for all the neighbours of whichever
    node the algorithm is considering at that point in time and hence we can safely assume that
    it will be of orfer approximatly close to the average degree number , however in the worst case
    where we have that all nodes are connected to each other directly then the order will be N again.
    We can now say that in a worst case scenario our algorithm will have a leading order of O(N^2),
    while it is also safe to assume that if we denote M as being the average degree then,
    the leading order in that case will be O(N*M).
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
        #print("***x=",x,' ***')
        for i in range(len(A[x])):              # Loop through the adjacency matrix
            v = A[x][i][0]                      # Get the first element of the tuples(node)
            if L2[v]==0 and a0*A[x][i][1] >= amin: # Condition from signal traffic
                Q.append(v)                     # Add unexplored neighbors to back of queue
                L2[v]=1
                L4[v].extend(L4[x])             # Add path to node x and node v to path
                L4[v].append(v)
            #print("v=",v)
            #print("Q=",Q)
    L5 = L4[J2]                                 # Printing the feasable path found

    return L5

#------------------------------------------------------------------------------#
#--------------------------------Part 1.2.2------------------------------------#
#------------------------------------------------------------------------------#

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

    Discussion: I will first discuss the implementation and the idea behind the
    algorithm before moving on to analyzing the time complexity of the code itself.

    a) Implementation: I based the implementation of my code largely on Dijkstra's algorithm
    which was covered in lectures. As we know this algorithm provides us with the shortest
    path in a network and since that is not what we are after we had to make some changes to
    the algorithm itself. First it is worth to note that while the original Dijakstra's
    implementation relies on weighted networks where the weight represents the distance and
    our goal is to minimize this weight or distance. We now find ourselves in the opposite
    situation, since we are talking about signal loss and the weights between nodes actually
    represent the min loss ratio that is admited, we will want to maximize this and hence
    our changes to the initial algorithm. Our first change is the switching of sign when
    checking if the weight is bigger than the initial weights which are initialized to a negative number
    in the same way we had previously distances to be initialized to "infinity".
    What we do now if fairly simple as we have taken the maximum value from the initial loop
    we will loop through the adjacency of the node with the max value and use this to find the
    maximum of minimums as explained before. Once this is done, we have only one issue at hand, which
    is dealing with the instance where the weight is 0 and where we get division by zero,
    this can easily be avoided with an if statement as its done. It is worth to note that
    the algorithm does not take in graphs with isolated nodes as the question clearly postulated
    that we are dealing with a problem where we have N junctions and hence no nodes will
    be possibly left disconnected.

    b) Time complexity: As per most of the algorithm thus far we start by initialization
    of the variables we require. We do so by initializing a dictionary which will contain the
    unexplored nodes, hence such operation will have O(N) time complexity where N is the len(A).
    We then immediately introduce a while loop that will run until there are no unexplored
    nodes left to explore. Assuming the worst scenario this will also be O(N). Inside such
    while loop we have two non-nested for loops , one of which will run for a worst case
    estimation of O(N) while the other which loops through the neighbours of a node, we can estimate
    will be running in the following time complexity O(M) where M is denoted to be the average
    degree number in our network.
    What we have then in terms of time complexity is the following:
    N+N*(N+M) = N + N^2 + N*M. Hence this will depend largely on the average number of nodes
    of our graph. however we may assume that when dealing with large graphs we will have
    N>>M which then leaves us with a leading order term of O(N^2) and hence our code
    has exponential time complexity. 
    """

    #----------------------Initialize dictionaries-----------------------------#

    a_init = -1                                 # Initial a value
    Edict = {}                                  # Explored nodes dict
    Udict = {}                                  # Unexplored nodes dict

    for n in range(len(A)):
        Udict[n] = a_init                       # Initializing a as -1
    Udict[J1] = 1                               # a needed to source is 1

    #-----------------------------Main Search----------------------------------#

    while len(Udict) > 0:                       # While there are unexplored nodes
        #Find node with max a in Udict and move to Edict
        a_max  = a_init                         # Setting maximum a to intial one for starting value
        for n,w in Udict.items():               # Looping through the unexplored dictionary
            if w > a_max:                       # If the weight node is greater than current a_max
                a_max = w                       # Update a_max
                n_max = n                       # Keep track of node where this occurs
        Edict[n_max] = Udict.pop(n_max)         # Mark such node as visited
        #print("moved node", n_max)

        # Update provisional a's for unexplored neighbours of n_max
        for i in range(len(A[n_max])):          # Looping through nodes and their associated weights
            n = A[n_max][i][0]                  # Defining node
            w = A[n_max][i][1]                  # Defining associated value a of the node
            if n in Udict:                      # While it is in the unexplored dictionary
                Udict[n] = max(min(w,a_max),Udict[n]) #update value to maximum of all values on path

    if Edict[J2] != 0 :                         #Avoid division by zero
        a0_min = amin/Edict[J2]
        path = findPath(A,a0_min,amin,J1,J2)
        output = a0_min, path
        if len(path) == 0:
            output = -1,[]
    else:
        output = -1,[]

    return output,Edict

#-----------------------------Main not used------------------------------------#
#if __name__=='__main__':
    #add code here if/as desired
    #L=None #modify as needed
