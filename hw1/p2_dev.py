"""M345SC Homework 1, part 2
Anas Lasri and CID:01209387
"""

def bsearch(L,x):
    #Set initial start and end indices for full list
    istart = 0
    iend = len(L) - 1

    #Iterate and contract "active" portion of list
    while istart <= iend:
        imid = int(0.5*(istart+iend))
        if x==L[imid]:
            return imid
        elif x < L[imid]:
            iend = imid-1
        else:
            istart = imid + 1

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

    #L is the list that has N sublists, So we start our search process by
    #considering it first. Since we are using numpy we will have to loop
    #through the sublists of L, *N*, one by one.
    A = []
    for idx,sublist in enumerate(L):
        if L[idx][P] < target:
            for i in range(P,len(L[idx])):
                if L[idx][i] == target:
                    A.append([idx,i])
        elif L[idx][P] >= target:
            for i in range(len(L[idx])):
                if L[idx][i] == target:
                    A.append([idx,i])



    return A


def nsearch_time():
    """Analyze the running time of nsearch.
    Add input/output as needed, add a call to this function below to generate
    the figures you are submitting with your codes.

    Discussion: (add your discussion here)
    """



    return None #Modify as needed


if __name__ == '__main__':

    #add call(s) to nsearch here which generate the figure(s)
    #you are submitting
    nsearch_time() #modify as needed
