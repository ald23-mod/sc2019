"""M345SC Homework 1, part 1
Anas Lasri Doukkali and CID:01209387
"""

def common_elements(list1,list2):
    result = []
    for i in list1:
        mutation_count = 0
        for j in list2:
            if j == i:
                mutation_count += 1
        result.append(mutation_count)
    return result

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
#------------------------------Part1.1--------------------------------------

    location = {}                          #location
    test = {}
    k_mers = {}                            #List containing result above
    A = {}
    final = {}

    for i in range(len(S)-k+1):       #Looping through the initial string
        k_mer = S[i:i+k]              #defining the k_mer
        if k_mer not in k_mers:
            location[k_mer] = i       #Gives the location of each k_mer
            test[i] = k_mer
            k_mers[k_mer] = 1
            A[k_mer] = [i]
        else:
            location[k_mer] = i
            k_mers[k_mer] += 1
            A[k_mer].append(location.get(k_mer))

    for i in range(len(S)-k+1):
        k_mer = S[i:i+k]
        if len(A[k_mer]) > f-1:
            final[k_mer] = A[k_mer]


#------------------------------Part1.2-------------------------------------

    #The idea is as follows: I will now go through the dictionary of k-mers
    #remove the xth component from the k-mer, making the element a (k-1)-mer
    #and now we can deploy the same method used in part 1.1
    freq_dict = []              #ammendable list of frequent k-mers
    all_kmers = []              #ammendable list of all k-mers

    for key, value in final.items():
        temp = key[:x] + key[x+1:]
        freq_dict.append(temp)

    for key,value in k_mers.items():
        temp = key[:x] + key[x+1:]
        all_kmers.append(temp)


    L1,L2,L3=[final.keys()],[final.values()],[common_elements(freq_dict,all_kmers)]

    return L1,L2,L3,k_mers





if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    S='CCCTATGTTTGTGTGCACGGACACGTGCTAAAGCCAAATATTTTGCTAGT'
    #k=3
    #x=2
    #f=2
    #L1,L2,L3=ksearch(S,k,f,x):
