import numpy as np


""" U = Id; V = Id; S = BD;
For i from 0 to NMax
    (Q1, R1) = decomp_qr(matrix_transpose(S))
    (Q2, R2) = decomp_qr(matrix_transpose(R1))
    S = R2;
    U = U * Q2;
    V = matrix_transpose(Q1) * V
    End for
Return (U,S,V) """

def QR_transfo(S,N):
    n,m = np.shape(S)
    U = np.eye(n)
    V = np.eye(n)
    for i in range(N):
        Q1, R1 = np.linalg.qr(np.transpose(S))
        Q2, R2 = np.linalg.qr(np.transpose(R1))
        S = R2
        U = np.dot(U,Q2)
        V = np.dot(np.transpose(Q1),V)
    return (U,S,V)

def QR_simple(S,N):
    for i in range (N):
        Q, R = np.linalg.qr(S)  
        S = np.dot(R,Q)
    return S


# Création d'une matrice bidiagonale de taille (3,4)
mat = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        if i == j:
            mat[i][j] = 2
        elif i == j-1 :
            mat[i][j] = 1

print(mat)
W=np.dot(QR_transfo(mat,100)[0],QR_transfo(mat,100)[1])
print(np.dot(W,QR_transfo(mat,100)[2]))
print("\n")
print(QR_transfo(mat,5)[1])
print(QR_simple(mat,100))
