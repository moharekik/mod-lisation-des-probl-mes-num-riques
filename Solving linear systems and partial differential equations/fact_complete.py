import numpy as np

"""
cholesky(A) : performs the cholesky factorization
"""

def cholesky(A):
    n = len(A)
    L = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1):
            s = 0
            for k in range(j):
                s += L[i,k] * L[j,k]
            if i == j:
                L[i,j] = np.sqrt(A[i,i] - s)
            else:
                L[i,j] = ((A[i,j] - s)/L[j,j])
    return L
