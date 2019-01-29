"""M345SC Homework 1, part 1
Anas Lasri Doukkali and CID:01209387
"""

def ksearch(S,k,f,x):
    """
    Search for frequently-occurring k-mers within a DNA sequence
    and find the number of point-x mutations for each frequently
    occurring k-mer.
    Input:
    S: A string consisting of A,C,G, and Ts
    k: The size of k-mer to search for (an integer)
    f: frequency parameter -- the search should identify k-mers which
    occur at least f times in S
    x: the location within each frequently-occurring k-mer where
    point-x mutations will differ.

    Output:
    L1: list containing the strings corresponding to the frequently-occurring k-mers
    L2: list containing the locations of the frequent k-mers
    L3: list containing the number of point-x mutations for each k-mer in L1.

    Discussion: Add analysis here
    """
    #I now need to define code that will find all k(function input parameter)
    #-mers that occur at least f(function input parameter) times.
    #The return of this first calculation should be the location of these 
    #k-mers as well as the k-mers themselves.

    k_mers = {}                            #List containing result above

    for i in range(len(S)-k+1):       #Looping through the initial string
        k_mer = S[i:i+k]
        if k_mer not in k_mers:
            k_mers[k_mer] = 1
        else:
            k_mers[k_mer] += 1
    return k_mers

    




    L1,L2,L3=[],[],[]

    return L1,L2,L3, count_mer




if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    S='CCCTATGTTTGTGTGCACGGACACGTGCTAAAGCCAAATATTTTGCTAGT'
    #k=3
    #x=2
    #f=2
    #L1,L2,L3=ksearch(S,k,f,x):
