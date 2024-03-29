import numpy as np
from math import *

def rectangle_left(f,a,b,n):
    h=(b-a)/n
    s=0
    for i in range(n):
        s+=f(a+i*h)
    return h*s

def rectangle_right(f,a,b,n):
    h=(b-a)/n
    return rectangle_left(f,a,b,n)+h*(f(a+n*h)-f(a))

def rectangle_middle(f,a,b,n):
    h=(b-a)/n
    s=0
    for i in range(n):
        s+=f(a+i*h+h/2)
    return h*s

def trapeze(f,a,b,n):
    h=(b-a)/n
    return rectangle_left(f,a,b,n)+(h/2)*(f(b)-f(a))

def simpson(f,a,b,n):
    h=(b-a)/n
    s1=0
    s2=f(a+h/2)
    for i in range(1,n):
        s1+=f(a+i*h)
        s2+=f(a+i*h+h/2)
    I=h*((f(a)+f(b))/6+s1*1/3+s2*2/3)
    return [s1,s2,I]

##doublement de pas pour les methodes rectangle_left rectangle_right trapeze

def double(method,f,a,b,n):
    h=(b-a)/n
    I=method(f,a,b,n)/2
    s=0
    for i in range(n):
        s+=f(a+h*i+h/2)
    return I+(h/2)*s

#triplement de pas pour la methode rectangle_middle

def triple(f,a,b,n):
    h=(b-a)/n
    I=rectangle_middle(f,a,b,n)
    s=0
    for i in range(n):
        s+=f(a+i*h+h/6)+f(a+i*h+5*h/6)
    return (I/3)+(h/3)*s


##doublement de pas pour la methode simpson

def double_simpson(f,a,b,n):
    h=(b-a)/n
    I0=(f(a)+f(b))/6
    L=simpson(f,a,b,n)
    s1=L[0]
    s2=L[1]
    I1=s1+s2
    I2=0
    for i in range(2*n):
        I2+=f(a+i*h/2+h/4)
    return (h/2)*(I0+I1/3+2*I2/3)


##Fonction d'integration generique qui prend en parametre la methode d'integration
##Retourne la valeur de l'integral + le nombre d'itérations

def integration(method,f,a,b,n,eps):
    if (method==simpson):
        I=simpson(f,a,b,n)[2]
        I2n=double_simpson(f,a,b,n)
        while (abs(I2n-I)>eps):
            n=2*n
            I=I2n
            I2n=double_simpson(f,a,b,n)
    elif (method==rectangle_middle):
        I=method(f,a,b,n)
        I2n=triple(f,a,b,n)
        while (abs(I2n-I)>eps):
            n=3*n
            I=I2n
            I2n=triple(f,a,b,n)      
    else:
        I=method(f,a,b,n)
        I2n=double(method,f,a,b,n)
        while (abs(I2n-I)>eps):
            n=2*n
            I=I2n
            I2n=double(method,f,a,b,n)
    return [I2n,n]

#test d'efficacite

def g(x):
    return x**2+1

print("f(x)=x^2+1")
print("rectangle_left",integration(rectangle_left,g,0,1,1,10e-4))
print("rectangle_right",integration(rectangle_right,g,0,1,1,10e-4))
print("rectangle_middle",integration(rectangle_middle,g,0,1,1,10e-4))
print("trapeze",integration(trapeze,g,0,1,1,10e-4))
print("simpson",integration(simpson,g,0,1,1,10e-4))

def g2(x):
    return sinh(x)

print("f(x)=sinh(x)")
print("rectangle_left",integration(rectangle_left,g2,0,1,1,10e-4))
print("rectangle_right",integration(rectangle_right,g2,0,1,1,10e-4))
print("rectangle_middle",integration(rectangle_middle,g2,0,1,1,10e-4))
print("trapeze",integration(trapeze,g2,0,1,1,10e-4))
print("simpson",integration(simpson,g2,0,1,1,10e-4))


def g3(x):
    return atan(x)

print("f(x)=atan(x)")
print("rectangle_left",integration(rectangle_left,g3,0,1,1,10e-4))
print("rectangle_right",integration(rectangle_right,g3,0,1,1,10e-4))
print("rectangle_middle",integration(rectangle_middle,g3,0,1,1,10e-4))
print("trapeze",integration(trapeze,g3,0,1,1,10e-4))
print("simpson",integration(simpson,g3,0,1,1,10e-4))

