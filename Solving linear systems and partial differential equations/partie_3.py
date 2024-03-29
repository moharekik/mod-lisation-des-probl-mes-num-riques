import numpy as np 
import matplotlib.pyplot as plt
from exercice2_1 import *



def space_discretization(N):
    x = np.linspace(0, 1,N)
    y = np.linspace(0,1,N)
    h = 1/(N-1)
    return (x,y,h)

def fo_partial_derivative(T,h,i,j) :
    return ( (T[i][j+1] - T[i][j-1])/2*h , (T[i+1][j] -T[i-1][j])/2*h )

def so_partial_derivative(T,h,i,j) : 
    return ( (T[i][j+1] - T[i][j-1] - 2*T[i][j])/2*h , (T[i+1][j] -T[i-1][j] - 2*T[i][j])/2*h )


def laplacien(N):
    x = space_discretization(5)
    Mc=[[0]*N**2 for i in range (N**2)]
    for i in range (N**2):
        Mc[i][i]=-4*(1/x[2]**2)
        if i<N**2-N :
            Mc[i][i+N]=1/x[2]**2
            Mc[i+N][i]=1/x[2]**2
        if i<N**2-1 :
            Mc[i][i+1]=1/x[2]**2
            Mc[i+1][i]=1/x[2]**2
    for j in range (1,N):
        Mc[N*j-1][N*j]=0
        Mc[N*j][N*j-1]=0
    for i in range (N**2):
        print(Mc[:][i])
    return Mc
 
N=20
def central (N):
    F =  np.zeros((N,N))
    F[N//2][N//2]=10
    B = np.reshape(F,(N**2,1))
    return B

def edge(N):
    F=np.zeros((N,N))
    for i in range (N):
        F[0][i]=10
    B = np.reshape(F,(N**2,1))
    return B

def diagonal(N):
    F=np.zeros((N,N))
    for i in range (N):
        F[i][i]=10
    B = np.reshape(F,(N**2,1))
    return B



def solve_eq(dist):
    a=laplacien(N)
    b=dist(N)
    x = np.zeros((N**2,1))
    T = gradient_conj(a, b, x)
    T=np.reshape(T,(N,N))
    return T


def solve_eq_predefined(dist):
    a=laplacien(N)
    b=dist(N)
    T = np.linalg.solve(a,b)
    T=np.reshape(T,(N,N))
    return T

def erreur_relative(dist):
    T1=solve_eq(dist)
    T2=solve_eq_predefined(dist)
    x=np.linalg.norm(T1-T2)/np.linalg.norm(T2)
    return x

e1=erreur_relative(central)
e2=erreur_relative(edge)
e3=erreur_relative(diagonal)

print("pour N =",N)
print(" l'erreur relative pour la distribution centrale égale à",e1)
print(" l'erreur relative pour la distribution latérale égale à",e2)
print(" l'erreur relative pour la distribution diagonale égale à",e3)

##fig = plt.figure(figsize=(N,N))
##plt.subplot(2,3,1)
##plt.imshow(solve_eq(central), cmap= plt.cm.hot,interpolation = 'bicubic')
##plt.title("distribution point central")
##plt.subplot(2,3,4)
##plt.imshow(solve_eq_predefined(central), cmap= plt.cm.hot,interpolation = 'bicubic')
##plt.subplot(2,3,2)
##plt.imshow(solve_eq(edge), cmap= plt.cm.hot,interpolation = 'bicubic')
##plt.title("distribution barre latérale")
##plt.subplot(2,3,5)
##plt.imshow(solve_eq_predefined(edge), cmap= plt.cm.hot,interpolation = 'bicubic')
##plt.subplot(2,3,3)
##plt.imshow(solve_eq(diagonal), cmap= plt.cm.hot,interpolation = 'bicubic')
##plt.title("distribution barre diagonale")
##plt.subplot(2,3,6)
##plt.imshow(solve_eq_predefined(diagonal), cmap= plt.cm.hot,interpolation = 'bicubic')
##plt.show()