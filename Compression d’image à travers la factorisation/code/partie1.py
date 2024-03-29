
import numpy as np
import time
import matplotlib.pyplot as plt

#---------------Question 1---------------------

def householder_vector(U, V):
    if np.all(U == V ):
        return np.zeros(U.shape)
    W = U - V
    norm_W = np.linalg.norm(W)
    return W / norm_W


def householder(U, V) :
    if np.allclose(U,V):
        return np.identity(np.shape(U)[0])
    N = householder_vector(U, V)
    H = np.identity(len(U)) - 2 * np.outer(N, N)
    return H


#--------------------- Test 1 ----
U=np.array([[3],[4],[0]])
V=np.array([[0],[0],[5]])

H = householder(U, V)
if __name__ == "__main__":
    print(H.dot(U)) #ca donne le vecteur V


#-------------- Question 2_1-------

def HouseHolderMulVector(U,V,x):
    N = householder_vector(U, V)
    opti = np.dot(N.T,x)
    return x-2*N*opti # le produit de n.T et x donne un scalaire, ce qui rend les choses simples.

#------------ Test 2------------


x = np.array([[1],[0], [0]])
H=householder(U,V)
exp_res =np.dot(H,x)
res = HouseHolderMulVector(U,V,x)
assert np.allclose(res, exp_res), f"Expected {exp_res}, but got {res}."

#---------------- Question 2_2--------

def HouseHolderMulMatrix(U,V,M):
    N = householder_vector(U, V)
    L = len(N)
    C = len(M[0])
    res = np.zeros((L,C))
    
    for i in range (C):
        nested_list = M[:, i:i+1]
        res[:,i] = np.ravel(HouseHolderMulVector(U, V, nested_list)) 
    return res
#----------------- Test 3---------------

 # Définir la matrice M à multiplier par H
M = np.array([[1, 1, 0], [1, 1, 0], [1, 0, 1]])
 
 # Calculer le produit HM
HM =  HouseHolderMulMatrix(U,V,M)
 
exp_HM = np.dot(H, M)

assert np.allclose(HM, exp_HM)

#----------- Graphe de comparaison de la complexite------

def mat_mul(A,B):
    C = np.zeros_like(A)
    n = len(A)
    for i in range (n):
         for j in range (n):
              for k in range (n):
                  C[i][j] += A[i][k]*B[k][j]
    return C

def random_matrix(n):
    return np.random.rand(n, n)

n_range = range(10, 101, 10)

time_dot = []
time_hh = []

for n in n_range:

    M = random_matrix(n)
    U = np.random.rand(n, 1)
    V = np.random.rand(n, 1)
    H = householder(U, V)
    
    start = time.time()
    mat_mul(H,M)
    end = time.time()
    time_dot.append(end - start)

    
    start = time.time()
    HouseHolderMulMatrix(U, V, M)
    end = time.time()
    time_hh.append(end - start)

plt.plot(n_range, time_dot, label="Dot Method")
plt.plot(n_range, time_hh, label="Householder Method")
plt.xlabel("Matrix Size")
plt.ylabel("Run Time (seconds)")
plt.legend()
plt.show()
