import math
import numpy as np
from partie1 import *

""" Fonction qui prend une matrice de HouseHolder et une taille de matrice pour créer la matrice Q en complétant si nécessaire avec la matrice identité. La matrice
de HouseHolder doit nécessairement être carrée. """
def hhPadding(hh_matrix, n):
    diff = n - np.shape(hh_matrix)[0]
    if diff < 0:
        return -1
    if diff == 0:
        return hh_matrix


    M = np.zeros((np.shape(hh_matrix)[0], diff))
    res = np.concatenate((M, hh_matrix), axis=1)

    id_matrix = np.identity(diff)
    N = np.zeros((diff, np.shape(hh_matrix)[0]))
    N = np.concatenate((id_matrix, N), axis=1)

    res = np.concatenate((N, res))
    return res

#Test de la fonction hhPadding
if __name__ == "__main__":
    n = 6
    M = np.array([[1, 2, 5],
        [3, 4, 7],
        [6, 9, 1]])
    print("Matrice de HouseHolder :\n", M)

    res = hhPadding(M, n)
    print("Résultat du remplissage avec une matrice finale de taille 6 :\n", res)
#Fin Test



def bidiagonale (BD) :
    Qleft = np.identity(np.shape(BD)[0])
    Qright = np.identity(np.shape(BD)[1]) 
    n=min(np.shape(BD)[0],np.shape(BD)[1])
    for i in range(n):
        # Let Q1 be a HH matrix mapping BD[i:n,i] on a vector with a single non-zero element
        U1=BD[i:np.shape(BD)[0],i]
        V1=np.zeros(np.shape(BD)[0]-i)
        Q1=householder(U1,V1)

        Q1 = hhPadding(Q1, n)

        ##la fonction qui cree la matrice HH est ID
        Qleft = np.dot(Qleft,Q1)
        #BD = HouseHolderMulMatrix(U1,V1,BD)
        BD = np.dot(Q1,BD)
        if (i < (np.shape(BD)[1] - 2)):
        # Let Q2 be a HH matrix mapping BD[i,(i+1):m] on a vector with a single non-zero element
            U2=BD[i+1:np.shape(BD)[0],i]
            V2=np.zeros(np.shape(BD)[0]-i-1)
            Q2 = householder(U2, V2)
            Q2 = hhPadding(Q2, n)
            Qright = np.dot(Q2,Qright)
            #BD = HouseHolderMulMatrix(U2,V2,np.transpose(BD))
            #BD=np.transpose(BD)
            BD=np.dot(BD,Q2)
    return (Qleft, BD, Qright)



if __name__ == "__main__":
    M = np.array([[2, 4, 2],[4, 1, 3],[1, 3, 2]])
    print("Matrice test de:\n", M)

    QL,BD,QR=bidiagonale(M)
    print(BD)
