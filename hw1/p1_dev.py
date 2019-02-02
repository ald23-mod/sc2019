"""M345SC Homework 1, part 1
Anas Lasri Doukkali and CID:01209387
"""
import numpy as np
import time
import matplotlib.pyplot as plt
nucleo = ['A','T','C','G']
def generate(N):
    seq = ''.join(list(np.random.choice(nucleo,N)))
    return seq

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
    k_mers = {}                            #List containing result above
    A = {}
    final = {}
    test = []

    for i in range(len(S)-k+1):       #Looping through the initial string
        k_mer = S[i:i+k]
        test.append(k_mer)              #defining the k_mer
        if k_mer not in k_mers:
            location[k_mer] = i       #Gives the location of each k_mer
            k_mers[k_mer] = 1
            A[k_mer] = [i]
            if len(A[k_mer]) > f-1:
                final[k_mer] = A[k_mer]

        else:
            location[k_mer] = i
            k_mers[k_mer] += 1
            A[k_mer].append(location.get(k_mer))
            if len(A[k_mer]) > f-1:
                final[k_mer] = A[k_mer]

#------------------------------Part1.2-------------------------------------

    #The idea is as follows: I will now go through the dictionary of k-mers
    #remove the xth component from the k-mer, making the element a (k-1)-mer
    #and now we can deploy the same method used in part 1.1
    all_kmers = []              #ammendable list of all k-mers
    freq_dict = []


    for key, value in final.items():
        temp = key[:x] + key[x+1:]
        freq_dict.append(temp)


    for str in test:
        temp = str[:x] + str[x+1:]
        all_kmers.append(temp)

    count  = []
    for i in range(len(freq_dict)):
        count.append(all_kmers.count(freq_dict[i]) - len(list(final.values())[i]))

    L1,L2,L3 = list(final.keys()),list(final.values()),count
    return L1,L2,L3

#-------------------------Part1.3/Analyze-----------------------------------
#In this part I will produce plots.
def analyze(k,f,x):
    final_dt = []
    size = []
    for  i in range(100000,1000000,100000):
        dt = []
        size.append(i)
        for j in range(10):
            t1 = time.time()
            dummy = ksearch(generate(i),k,f,x)
            t2 = time.time()
            dt.append(t2-t1)
        final_dt.append(sum(dt)/len(dt))

    plt.plot(size[:],final_dt[:],'x--')

if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    S='CCCTATGTTTGTGTGCACGGACACGTGCTAAAGCCAAATATTTTGCTAGT'
    #k=3
    #x=2
    #f=2
    #L1,L2,L3=ksearch(S,k,f,x):
