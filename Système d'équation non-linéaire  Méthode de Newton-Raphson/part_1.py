import numpy as np

def Newton_Raphson(f,J,U0,N,epsilon):
    Xn=U0
    if (type(U0)==float)or (type(U0)==int):
        for i in range(N):
            Xn=Xn-f(Xn)/J(Xn)
            if abs(Xn)<epsilon:
                return Xn
        return Xn

    for i in range(N):
        h=np.linalg.lstsq(J(Xn),f(Xn),rcond=-1)
        Xn=np.subtract(Xn,h[0])
        if np.linalg.norm(Xn)<epsilon:
            return Xn
    return Xn



'''QUESTION 4'''


def Newton_Raphson_backtrack(f,J,U0,N,epsilon):

    Xn=U0
    if (type(U0)==float)or (type(U0)==int):
        for i in range(N):
            h=f(Xn)/J(Xn)
            Xn_1=Xn-h
            if f(Xn_1)-f(Xn)>0:
                Xn_1=Xn-0.5*h
            else:
                Xn_1=Xn-2*h


            if abs(Xn_1)<epsilon:
                return Xn
            Xn=Xn_1
        return Xn



    for i in range(N):


        h=np.linalg.lstsq(J(Xn),f(Xn),rcond=-1)
        Xn_1=np.subtract(Xn,h[0])
        if np.linalg.norm(np.subtract(f(Xn_1),f(Xn)))>0:
            h[0][1]=2*h[0][1]
            h[0][0]=2*h[0][0]
            Xn_1=np.subtract(Xn,h[0])
        else:
            h[0][1]=0.5*h[0][1]
            h[0][0]=0.5*h[0][0]
            Xn_1=np.substract(Xn,h[0])
        if np.linalg.norm(Xn)<epsilon:
            return Xn
        Xn=Xn_1
    return Xn

def poly(x):
    return x**3-x**2

def deri_poly(x):
    return 3*x**2-2*x

