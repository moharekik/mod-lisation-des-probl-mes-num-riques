import numpy as np
import random 
import fact_complete as fc



def def_po(A):
    N=A.shape[0]
    for i in range(N):
        sum=0
        for j in range(N):
            sum+=A[i][j]
        A[i][i]=sum+1


def mat_sym_creuse(n,p):#matrix n*n, p zeros. p must be even
    A=np.zeros([n,n])
    compteur =0
    coord=[]
    A[0][0]=random.randint(1,101)
    A[n-1][n-1]=random.randint(1,101)
    for k in range(1,n):#add numbers, symetrically
        for i in range(k):
           randintt=random.randint(1,101)
           A[i][k]=randintt
           A[k][i]=randintt

    for i in range(p//2):#add zeros
        randintt=random.randint(0,n-1)
        randintt2=random.randint(0,n-1)#random coords
        lsite=[randintt,randintt2]
        lsite2=[randintt2,randintt]
        while ((lsite in coord) or (lsite2 in coord) or (randintt==randintt2)):#dont add zeros where there already are zeros
            randintt=random.randint(0,n-1)
            randintt2=random.randint(0,n-1)
            lsite=[randintt,randintt2]
            lsite2=[randintt2,randintt]

        coord.append([randintt,randintt2])
        coord.append([randintt2,randintt])
        A[randintt][randintt2]=0
        A[randintt2][randintt]=0
    def_po(A)
    return A



"""
cholesky(A) : performs the incomplete cholesky factorization
"""

def cholesky_incomplet(A):
    n = len(A)
    L = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1):
            s = 0
            if A[i][j]==0:
                L[i][j]=0.
            else :
                for k in range(j):
                    s += L[i,k] * L[j,k]
                if i == j:
                    L[i,j] = np.sqrt(A[i,i] - s)
                else:
                    L[i,j] = ((A[i,j] - s)/L[j,j])
    return L        




