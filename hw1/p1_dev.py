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

    Discussion:

    """

#------------------------------Part1.1--------------------------------------

#Initialization of Data Structures:------------------------

    k_mers = {}                            #Dictionary containing the k-mers
    A = {}                                 #Dictionary containing loc. for repeats
    final = {}                             #Dictionary the final result, 1.1

    size = len(S)-k+1                      #Precalculating the size for the loops

#-----------------------------------------------------------

    for i in range(size):                 #Looping through the initial string
        k_mer = S[i:i+k]                  #defining the k_mer at each iteration
        if k_mer not in k_mers:
            k_mers[k_mer] = 1
            A[k_mer] = [i]
            if len(A[k_mer]) >= f:
                final[k_mer] = A[k_mer]
        else:
            k_mers[k_mer] += 1
            A[k_mer].append(i)
            if len(A[k_mer]) >= f:
                final[k_mer] = A[k_mer]

    L1 = list(final.keys())
    L2 = list(final.values())
#------------------------------Part1.2-------------------------------------

    #The idea is as follows: I will now go through the dictionary of k-mers
    #remove the xth component from the k-mer, making the element a (k-1)-mer
    #and now we can deploy the same method used in part 1.1
    L3 = []
    var = {}

    for i in range(size):
        k_mer = S[i:i+x] + S[i+x+1:i+k]
        if k_mer not in var:
            var[k_mer] = 0
        if k_mer in var:
            var[k_mer] += 1
    for i in range(len(L1)):
        k_mer = L1[i][:x] + L1[i][x+1:]
        if k_mer in var:
            L3.append(var[k_mer] - len(L2[i]))

    return L1,L2,L3

#-------------------------Part1.3/Analyze-----------------------------------
#In this part I will produce plots.
#def analyze(k,f,x):
#    final_dt = []
#    size = []
#    for  i in range(2000,10000,100):
#        dt = []
#        size.append(i)
#        for j in range(10):
#            t1 = time.time()
#            dummy = ksearch(generate(i),k,f,x)
#            t2 = time.time()
#            dt.append(t2-t1)
#        final_dt.append(sum(dt)/len(dt))
#
#    plt.plot(size[:],final_dt[:],'x--')
#
#def analyze2(f,x):
#    final_dt = []
#    size = []
#    L = generate(100000)
#    for i in range(1,110000,5000):
#        dt = []
#        size.append(i)
#        t1 = time.time()
#        dummy = ksearch(L,i,f,x)
#        t2 = time.time()
#        dt.append(t2-t1)
#        final_dt.append(sum(dt)/len(dt))
#        plt.plot(size[:],final_dt[:],'x--')
#
#def analyze3(k,x):
#    final_dt = []
#    size = []
#    L = generate(100000)
#    for i in range(1,110000,5000):
#        dt = []
#        size.append(i)
#        t1 = time.time()
#        dummy = ksearch(L,k,i,x)
#        t2 = time.time()
#        dt.append(t2-t1)
#        final_dt.append(sum(dt)/len(dt))
#        plt.plot(size[:],final_dt[:],'x--')
#
#def analyze4(k,f):
#    final_dt = []
#    size = []
#    L = generate(100000)
#    for i in range(1,110000,5000):
#        dt = []
#        size.append(i)
#        t1 = time.time()
#        dummy = ksearch(L,k,f,i)
#        t2 = time.time()
#        dt.append(t2-t1)
#        final_dt.append(sum(dt)/len(dt))
#        plt.plot(size[:],final_dt[:],'x--')



if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    S='CCCTATGTTTGTGTGCACGGACACGTGCTAAAGCCAAATATTTTGCTAGT'
    #k=3
    #x=2
    #f=2
    #L1,L2,L3=ksearch(S,k,f,x):
