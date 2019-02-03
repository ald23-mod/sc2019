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
    3.1) Concise description of my algorithm:
    I start my algorithm by initializing the necessary dictionaries, these
    dictionaries are the following:
         -k_mers:This dictionary saves each k-mer as a key along with the values
         being the number of occurences of such k-mer. Thus, for example a
         k-mer that has occured 3 times will be in the dictionary under the
         form ('k-mer': 3)
         -A: This dictionary will have as the value of each k-mer all of its
         positions if it does appear more than once in a list
         -final: This is a dictionary with the most frequent k-mers along
         all the locations in which they can be found
    The algorithm is fairly simple and easy to understand as after I initialize
    the lists I will be using I start by looping through the string of DNA
    and considering every k-mer that arises by first assigning the k-mer to
    list named test,
    """

#------------------------------Part1.1--------------------------------------

#Initialization of Data Structures:------------------------

    k_mers = {}                            #Dictionary containing the k-mers
    A = {}                                 #Dictionary containing loc. for repeats
    final = {}                             #Dictionary the final result, 1.1
    test = []                              #List containing the k_mers

#-----------------------------------------------------------

    for i in range(len(S)-k+1):            #Looping through the initial string
        k_mer = S[i:i+k]                   #defining the k_mer
        test.append(k_mer)
        if k_mer not in k_mers:
            k_mers[k_mer] = 1
            A[k_mer] = [i]
            if len(A[k_mer]) > f-1:
                final[k_mer] = A[k_mer]
        else:
            k_mers[k_mer] += 1
            A[k_mer].append(i)
            if len(A[k_mer]) > f-1:
                final[k_mer] = A[k_mer]

#------------------------------Part1.2-------------------------------------

    #The idea is as follows: I will now go through the dictionary of k-mers
    #remove the xth component from the k-mer, making the element a (k-1)-mer
    #and now we can deploy the same method used in part 1.1
    all_kmers = []              #ammendable list of all k-mers
    freq_dict = []
    count  = []


    for key, value in final.items():
        temp = key[:x] + key[x+1:]
        freq_dict.append(temp)


    for str in test:
        temp = str[:x] + str[x+1:]
        all_kmers.append(temp)

    for i in range(len(freq_dict)):
        count.append(all_kmers.count(freq_dict[i]) - len(list(final.values())[i]))

    L1,L2,L3 = list(final.keys()),list(final.values()),count
    return L1,L2,L3,final

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

def analyze2(f,x):
    final_dt = []
    size = []
    L = generate(100000)
    for i in range(1,100000,5000):
        dt = []
        size.append(i)
        t1 = time.time()
        dummy = ksearch(L,i,f,x)
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
