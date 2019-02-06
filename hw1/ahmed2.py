
def bsearch(L,x):

    #Set initial start and end indices for full list
    istart = 0
    iend = len(L)-1

    #Iterate and contract "active" portion of list
    while istart<=iend:

        imid = int(0.5*(istart+iend))

        if x==L[imid]:
            while L[imid - 1] == L[imid]:
                imid = imid-1
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
    Lout=[]

    for i in range(len(L)):
        if L[i][P-1] >= target:
            j = bsearch(L[i][0:P],target)
            if j != -1000:
                Lout.append([i,j])
                while L[i][j] == L[i][j+1] and j+1 < P:
                    j = j + 1
                    Lout.append([i,j])

        for k in range(P,len(L[i])):
            if L[i][k] == target:
                Lout.append([i,k])

    return Lout
