import numpy as np
import fact_complete as fc
import time
import exercice1_2 as ifc
import matplotlib.pyplot as plt

A = np.array([[1, 1, 1, 1], [1, 5, 5, 5], [1, 5, 14, 14], [1, 5, 14, 15]])

B = np.array([[1,1,1],[1,2,2],[1,2,3]])

print("=========================")
print("TEST : cholesky()")
print("=========================\n")

print("-------------------------")
print("Matrix A : ")
print(A)
print("\n")

print("cholesky(A) : ")
print(fc.cholesky(A))
print("\n")

print("Vérification : A = T⋅T^t ")
print(np.dot(fc.cholesky(A),np.transpose(fc.cholesky(A))))
print("-------------------------")
print("\n\n")

print("-------------------------")
print("Matrix B : ")
print(B)
print("\n")

print("cholesky(B) : ")
print(fc.cholesky(B))
print("\n")

print("Vérification : B = T⋅T^t ")
print(np.dot(fc.cholesky(B),np.transpose(fc.cholesky(B))))
print("-------------------------")
print("\n\n")


print("=========================")
print("TEST : comparison between incomplete cholesky and complete cholesky")
print("=========================\n")

time_complete=[]
time_incomplete=[]
x=[]
dimension=30

n=(int(dimension*(dimension-1)))//2
print("There are 200 tests ,with a total of ", 2*n,"zeros in the matrix at the end")
for k in range(0,n,n//100):
    x.append(k)
    compteur_temps=0
    for i in range(20):
        A=ifc.mat_sym_creuse(dimension,2*k)
        start=time.time()
        ifc.cholesky_incomplet(A)
        end=time.time()
        timme=end-start
        compteur_temps+=timme
    time_incomplete.append(compteur_temps)
    if k%50==0:
        print("time for incomplete factorisation with ",2*k," zeros is ", timme,"\n")


for k in range(0,n,n//100):
    compteur_temps=0
    for i in range(20):
        A=ifc.mat_sym_creuse(dimension,2*k)
        start=time.time()
        fc.cholesky(A)
        end=time.time()
        timme=end-start
        compteur_temps+=timme
    time_complete.append(compteur_temps)
    if k%50==0:
        print("time for complete factorisation with ",2*k," zeros is ", timme,"\n")
    
plt.plot(x,time_complete,label="temps factorisation complete")
plt.plot(x,time_incomplete,label="temps factorisation incomplete")
plt.legend()
plt.show()



def conditionnement(A):
    T=ifc.cholesky_incomplet(A)
    Tprime=fc.cholesky(A)
    Tinv=np.linalg.inv(T)
    Tprimeinv=np.linalg.inv(Tprime)
    M=np.dot(np.transpose(Tinv),Tinv)
    Mprime=np.dot(np.transpose(Tprimeinv),Tprimeinv)
    cond=np.linalg.cond(np.dot(np.linalg.inv(M),A))
    condprime=np.linalg.cond(np.dot(np.linalg.inv(Mprime),A))
    return (cond,condprime)

print("chargement" )
compteur=[0,0,0]

for i in range(1000):
    A=A=ifc.mat_sym_creuse(10,60)
    if i%300==0:
        print(".")
    compteur[0]+=conditionnement(A)[0]#incomplet
    compteur[1]+=conditionnement(A)[1]#complet
    compteur[2]+=np.linalg.cond(A)

print("\n")

for k in range(3):
    compteur[k]=compteur[k]/1000


print("conditionnement de A : ", compteur[0], "conditionnement de T^(-1)(t)*T^(-1) pour cholesky complet : ",compteur[2])
print("conditionnement de T^(-1)(t)*T^(-1) cholesky incomplet : ",compteur[1])