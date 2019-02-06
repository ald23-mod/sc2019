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
    3.a) My algorithm is based completely in the dictionary data structure
    implemented in python. This data structure which was introduced in Lab2 is
    the most efficient method of all the different attempts made at this question.
    Another fruitful yet less efficient attempt at this question was the use of
    the Rolling Hash however my implementation of such algorithm turned out to be
    slower and it inevitably required the use of the inbuilt Python dictionaries
    which made my implementation redundant.

    Since dictionaries in Python are themselves implemented using Hash Tables and
    optimized I decided to implement Question 1 using mainly this data structure.
    The calculation for L1 and L2 was fairly simple and it required the initialization
    of two dictionaries which I called 'A' and 'final'. A stores as values all the
    k-mers and as keys to the values the position of such k-mers. final on the
    other hand is a subset of the first dictionary that is created to contain the
    frequent k-mers and their position.

    L3 calculation turned out to be trickier however I reached a optimized
    way using two separate for loops that reduces the time complexity. This is done
    by first using a loop to get the x-mutations of all the k-mers and puting that
    into a dictionary which I named var along with their count. Secondly I used
    another loop in which I made the final calculations which is basically substracting
    the value for each mutation key to the repetitions which we get from L2.

    3.b) To analyze the running time of my code I decided to do while going through
    the function and evaluating the time complexity at each loop while also using
    plots to confirm this. The code for the plots will be left in the .py file
    commented and hence the marker can decide whether or not to make use of it.
    Let us start our time complexity analysis I will now go through the entire
    algorithm pointing out the time order for each loop and operation in the hopes
    that at the end I will get a correct estimate for the timing.
    My code has a total of three for loops. The two first for loops are looping
    through the same length(len(S)-k+1=N).

    I will now analyze the operations inside this first loop:
    -The loop has size N, and hence I can say than it is in order O(N), but I still
    have to analyze the operations inside the loop themselves which will give us
    the leading term coefficient.
    -We start with an assigment of the k-mer to a dummy variable using slices (O(k),size of k-mer)
    called k-mer and we inmediately introduce an if/else statement. The if statement
    itself is O(1) and hence we will need to evaluate a worst/best case scenario.
    The best case scenario is when the if statement is true and hence the order inside
    is only that of assigning(O(1)). We won't have to worry about the nested if statement
    as this will only happen when f = 1 which is unlikely to happen. Hence the best
    case scenario for my first loop is O(2kN). The worst case happens if my loop follows
    the (else :) condition as it will then have to go through another assignment
    statement which makes my worst scenario of order O(3kN).

    I then move on to assign the values for L1 and L2 which were calculated during
    the previous loop. Each of these assignments is O(1), and hence of total, O(2).

    I now start the second loop which is also of size N and hence we will now need
    to analyze the operations inside said loop again:
    -This second loop starts with slicing which order O(final-start) in our case
    is O(x) + O(k-x-1) which is ~ O(k) and then proceeds with and if and else statement
    both with the same order as the same amount of work is being done by the code.
    Making the order time for this entire loop ~O((k+1)N).

    Moving on now to the last loop in our code. It is a for loop over the frequent
    k-mers which means its length will depend on how many frequent element we get,
    let us assume a worst case scenario in which this number tends to N. This along
    of everything we have said so far about slicing order makes the order of this entire loop
    O((3+k)N).

    All of the above calculations along with the fact that N = len(S)-k+1, leads
    us to the following conclusion. Let len(S) = n then the order of my entire function
    is approx. ~ (kn-k^2) where n is the len(S) and k is the length of the k-mers.
    All of this can be backed by the figures that I generated to confirm my line of
    thought.
    Finally I can claim that my function runs O(n)(linear), with leading order coeff.
    ~(4k-k^2/n) approximately. My analyze2 function plots a graph in a shape of a parabola
    that confirms the shape of the leading coefficient. 

    """

#------------------------------Part1.1--------------------------------------

#Initialization of Data Structures:-----------------------------------------

    A = {}                                 #Dictionary containing loc. for repeats
    final = {}                             #Dictionary the final result, 1.1

    size = len(S)-k+1                      #Precalculating the size for the loops

#Loop to define useful dictionaries----------------------------------------

    for i in range(size):                 #Looping through the initial string
        k_mer = S[i:i+k]                  #Defining the k_mer at each iteration
        if k_mer not in A:           #If statement set up
            A[k_mer] = [i]                #Adding position of the k-mer
            #The next line is redundant as length is one, however f may be set to 1
            if len(A[k_mer]) >= f:        #Checking for length
                final[k_mer] = A[k_mer]   #If statement holds then add to final
        else:                             #If k_mer is already present in dict.
            A[k_mer].append(i)            #Append location
            if len(A[k_mer]) >= f:        #Checking for length condition as above
                final[k_mer] = A[k_mer]   #Add to final result


    L1 = list(final.keys())               #Define L1
    L2 = list(final.values())             #Define L2
#------------------------------Part1.2-------------------------------------

#Initialization of Data Structures:----------------------------------------

    L3 = []                               #Initialize empty L3 list
    var = {}                              #Initialize dictionary

#Loop to define the dictionary from which L3 will be calculated------------

    for i in range(size):
        k_mer = S[i:i+x] + S[i+x+1:i+k]   #Defining the x-mutations from the k-mers
        if k_mer not in var:              #If statement to set counter
            var[k_mer] = 0
        if k_mer in var:                  #Setting count
            var[k_mer] += 1
    for i in range(len(L1)):              #Using L1 to loop through and calculate L3
        k_mer = L1[i][:x] + L1[i][x+1:]   #x-mutations of frequent k-mers
        if k_mer in var:
            L3.append(var[k_mer] - len(L2[i])) # Appending the result
    return L1,L2,L3

#-------------------------Part1.3/Analyze-----------------------------------
#In this part I will produce plots.
#def analyze(k,f,x):
#    final_dt = []
#    size = []
#    for  i in range(2000,50000,1000):
#        dt = []
#        size.append(i)
#        for j in range(10):
#            t1 = time.time()
#            dummy = ksearch(generate(i),k,f,x)
#            t2 = time.time()
#            dt.append(t2-t1)
#        final_dt.append(sum(dt)/len(dt))

#    plt.plot(size[:],final_dt[:],'x--')
#
#def analyze2(f,x):
#    final_dt = []
#    size = []
#    L = generate(10000)
#    for i in range(0,10000,800):
#
#        dt = []
#        size.append(i)
#        t1 = time.time()
#        dummy = ksearch(L,i,f,x)
#        t2 = time.time()
#        dt.append(t2-t1)
#        final_dt.append(sum(dt)/len(dt))
#        plt.plot(size[:],final_dt[:],'x--')

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
