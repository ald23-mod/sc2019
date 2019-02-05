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



    L1,L2,L3=[],[],[]
    X,Y,Z = {},{},[]

    for i in range(len(S) - k + 1):

        w = S[i:i+k]
        Z.append(w)
        if w not in X:
            X[w] = 1
            Y[w] = [i]

        else:
            X[w] = X[w]+1
            Y[w].append(i)

    for w in X:
        if X[w] >= f:
            L1.append(w)
            L2.append(Y[w])



    L = L1[:]

    #for r in range(len(L)):
        #L[r] = L[r][:x] + L[r][x+1:]

    #for z in range(len(Z)):
        #Z[z] = Z[z][:x] + Z[z][x+1:]


    for z in range(len(Z)):
        Z[z] = Z[z][:x] + Z[z][x+1:]
        if z in range(len(L)):
            L[z] = L[z][:x] + L[z][x+1:]



    for l in range(len(L)):
        L3.append(Z.count(L[l]) - len(L2[l]))




    return L1,L2,L3
