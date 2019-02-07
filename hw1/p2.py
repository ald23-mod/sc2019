"""M345SC Homework 1, part 2
Anas Lasri and CID:01209387
"""
import numpy as np
import time
import matplotlib.pyplot as plt


def generate2(N,M,P):
    """This is the function used to generate the randomized Initial
    L which is the main input into our nsearch function
    """
    A = np.random.randint(100000, size = (N,M))
    A[:,0:P] = np.sort(A[:,0:P])
    return A

#This is an edited version of the binary search in which it find the first
#ocurrence  of the target. This will help with the implementation later on
def bsearch(L,x):
    """This is the edited binary search as described below. I edited the function
    from the lectures to make it find the first instance of occurence of the
    target and so being able to add a while loop to find all occurences of it.
    """
    #Set initial start and end indices for full list
    istart = 0
    iend = len(L)-1
    #Iterate and contract "active" portion of list
    while istart<=iend:
        imid = int(0.5*(istart+iend))
        if x==L[imid]:
            while L[imid] == L[imid-1]: #Modified part to find the first ocurrence
                imid = imid - 1
            return imid
        elif x < L[imid]:
            iend = imid-1
        else:
            istart = imid+1

    return -1000

def nsearch(L,P,target):
    """Input:
    L: list containing *N* sub-lists of length M. Each sub-list
    contains M numbers (floats or ints), and the first P elements
    of each sub-list can be assumed to have been sorted in
    ascending order (assume that P<M). L[i][:p] contains the sorted elements in
    the i+1th sub-list of L
    P: The first P elements in each sub-list of L are assumed
    to be sorted in ascending order
    target: the number to be searched for in L

    Output:
    Lout: A list consisting of Q 2-element sub-lists where Q is the number of
    times target occurs in L. Each sub-list should contain 1) the index of
    the sublist of L where target was found and 2) the index within the sublist
    where the target was found. So, Lout = [[0,5],[0,6],[1,3]] indicates
    that the target can be found at L[0][5],L[0][6],L[1][3]. If target
    is not found in L, simply return an empty list (as in the code below)
    """

    A = []                                         #Initializing empty list
    for idx,sublist in enumerate(L):               #Looping over list
        if L[idx][P] < target:                     #Linear search
            for i in range(P,len(L[idx])):
                if L[idx][i] == target:           #Check whether its target
                    A.append([idx,i])             #Append to result list
        else:                                     #binary search
            j = bsearch(L[idx][0:P],target)
            if j != -1000:
                A.append([idx,j])
                while L[idx][j] == L[idx][j+1] and j+1 < P:  #Check for other target hits
                    j = j + 1
                    A.append([idx,j])
            for id in range(P,len(L[idx])):       #linear search
                if L[idx][id] == target:          #Check whether its target
                    A.append([idx,id])            #Append to result list


    return A


def nsearch_time():
    """Analyze the running time of nsearch.
    Add input/output as needed, add a call to this function below to generate
    the figures you are submitting with your codes.

    Discussion:
    2.a) I start my algorithm by looping thorugh the imput list L using
    enumerate which can return to us the index along with the sublist to deal with.
    I then inmediately introduce an if statement. This is done for the following reason:
    I will first go to the last sorted in each sublist of L and check wheter that element
    is less than the target, the reason for this is if said element is less than the
    target then we now for a fact that our target will not lie within the first P
    sorted elements and hence we can skip doing any search method on this part of
    our sub-lists. If this is the case I simply perform linear search on the undordered
    part of the sublist to find any matches with the target and add them to my result
    list.
    Let us now assume that the P'th element of the sublist is not less than the target,
    what my code will do then is treat the sublist as two separate problems, performing
    an edited version of binary search on the first half as it is ordered and then it
    will simply perform linear search in the second unordered part of the sublist.
    The last part to clarify about my algorithm is the binary search. The basic function
    was copied from the one used in the lectures with minor tweaks. the tweaks are to
    solve a minor issue which is the fact that binary search will only find one instance
    of the target in our sublist, and hence if the target is repeated within the sorted part
    it will not be able to identify the repeats. To solve this, the binary search
    function is edited to find the first instance or ocurrence of our target in the
    sorted sublist, and then from there within our main function we can add a while
    statement that will be able to count over all the ocurrences.

    2.b) As you can see from the algorithm we start with an enumerate loop of size N, where
    N is the number of sublists inside L itself. I then introduce an if statement as
    explained before. The best case scenario time-wise follows if the if statement is true, in
    which case we only use linear search on the unsorted part of the sublist, this
    as the name indicated is in linear order as we will need to loop through all the
    elements.
    Let us now consider the worst case for now which is when the algorithm has to
    do both binary and linear search. We know the order of binary search is logarithm
    of P in the base of 2 and since we already know from the previous part that linear
    search has linear order in (M-P) which is the length of the unsorted sublist.

    And hence we can say than in general the order of the function is:
    N*(log_{2}P + M-P) and hence  big-O is O(N(M-P)). As I discussed before  in the best case
    where we do not go through binary search we have order being N*(M-P). This is
    while counting the order without caring for leading constants.

    2.c) My algorithm can be considered as efficient because I first of all only consider
    the sorted array if the condition mentioned above holds and hence I do not even
    consider such part of the sublists which affects the time hugely. Secondly I
    even when I consider the sorted part of the list I created an edited version of the
    binary search that reduces the time and helps me to not employ linear search
    throughout all of the sublists

    2.e) I will now comment on my plots and relate them to the results I discussed
    theoretically. First I would like to say that for better time estimates I averaged
    over 10 simulations in all the plots. Let us discuss the plots one by one:

    fig1.png: To create this plot I changed the size of N in increasing order while
    keeping both M and P constant. Following the time complexity analysis that I completed
    this should have linear order as I said that the order nsearch follows the following
    expression N*(log_{2}P + M-P). This is exactly what the plot outputs.

    fig2.png: For this plot as you can see from the title of the plot I changed the M
    while keeping everything else constant and following the same theoretical principles
    from which we analysed the first plot this will also be linear in time and that
    is exactly what we get.

    fig3.png: For this plot I varied both N and M to find what would be the relation
    betqeen both parameters. This relation can be deduced from the expresion N*(log_{2}P + M-P)
    Since we will be varying M and N the important term will be N*(M-P) and since
    both of them are increasing in our plot we will get  a quadratic shape plot as the
    term N*M dictates.
    """

    timef,time1f, time2f = [],[],[]
    actual_time, actual_time1, actual_time2 = [],[],[]
    size, size1, size2 = [],[],[]

    for i in range(1000,10000,250):
        size.append(i)
        for j in range(10):
            L = generate2(i,100,50)
            t1 = time.time()
            nsearch(L,50,5)
            t2 = time.time()
            dt = t2 - t1
            timef.append(dt)
        actual_time.append(sum(timef)/len(timef))

    for i in range(1000,2000,100):
        N = 2000
        size1.append(i)
        for j in range(10):
            L = generate2(N,i,i//2)
            t1 = time.time()
            nsearch(L,i//2,50)
            t2 = time.time()
            dt = t2 - t1
            time1f.append(dt)
        actual_time1.append(sum(time1f)/len(time1f))

    for i in range(0,3000,250):
        size2.append(i)
        for j in range(10):
            L = generate2(i,i//2,i//4)
            t1 = time.time()
            nsearch(L,i//4,50)
            t2 = time.time()
            dt = t2 - t1
            time2f.append(dt)
        actual_time2.append(sum(time2f)/len(time2f))


    plt.figure()
    plt.plot(size[:],actual_time[:])
    plt.xlabel('Size of L')
    plt.ylabel('Time')
    plt.title('Anas Lasri, CID: 01209387 \n Plot varying the size of L')
    plt.savefig('fig1.png',dpi=400)

    plt.figure()
    plt.plot(size1[:],actual_time1[:])
    plt.xlabel('Size of M')
    plt.ylabel('Time')
    plt.title('Anas Lasri, CID: 01209387 \n Plot varying the size of M')
    plt.savefig('fig2.png',dpi=400)

    plt.figure()
    plt.plot(size2[:],actual_time2[:])
    plt.xlabel('Size of N')
    plt.ylabel('Time')
    plt.title('Anas Lasri, CID: 01209387 \n Plot varying the size of N and M')
    plt.savefig('fig3.png',dpi=400)


    return None #Modify as needed


if __name__ == '__main__':

    #add call(s) to nsearch here which generate the figure(s)
    #you are submitting
    nsearch_time() #modify as needed
