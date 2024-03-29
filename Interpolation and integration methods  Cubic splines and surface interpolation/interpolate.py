import sys
import numpy as np
from airfoil import load_foil
import matplotlib.pyplot as plt


# Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., y i = f ( x i ), with
# x 1 < x 2 < . . . < x N , and given values yp1 and ypn for the first derivative of the interpolating
# function at points 1 and n , respectively, this routine returns an array y2[1..n] that contains
# the second derivatives of the interpolating function at the tabulated points x i . If yp1 and/or
# ypn are equal to 1 × 10^30 or larger, the routine is signaled to set the correspond
def spline(xv,yv,yp1,ypn):
    n = len(xv)
    u = np.zeros(n-1)
    y2 = np.zeros(n)
    sig = 0
    p = 0
    if (yp1>1e30):
        y2[0] = 0       #the lower boundary condition is set either to be "natural"
        u[0] = 0        
    else:
        y2[0] = -1/2    #or else to have a specified first derivative
        u[0] = (3/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1)
    for i in range(1,n-1):                      #This is the decomposition loop of the tridiagonal algorithm.
        sig = (xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]) #y2 and u are used for temporary storage of the decomposed factors.
        p = sig*y2[i-1]+2
        y2[i] = (sig-1)/p
        u[i]=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1])
        u[i]=(6*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p
    if (ypn>1e30):
        qn = 0          #The upper boundary condition is set either to be "natural"
        un = 0
    else:               #or else to have a specified first derivative
        qn = 1/2
        un = (3/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]))
    y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1)
    k = n-2
    while (k>=0):
        y2[k]=y2[k]*y2[k+1]+u[k]    #this is the backsubstitution loop 
        k -= 1                      #of the tridiagonal algorithm
    return y2

# Given the arrays xa[1..n] and ya[1..n] , which tabulate a function (with the xa i ’s in order),
# and given the array y2a[1..n] , which is the output from spline above, and given a value of
# x , this routine returns a cubic-spline interpolated value y .


def splint(xa, ya, y2a, n, x):
    klo = 0
    khi = n-1
    while (khi-klo > 1):
        k = (khi+klo)//2
        if (xa[k] > x):
            khi = k
        else:
            klo = k
    h = xa[khi]-xa[klo]
    if (h == 0):
        print("Bad xa input to routine splint")
    a = (xa[khi]-x)/h
    b = (x-xa[klo])/h
    y1 = (ya[khi]-ya[klo])/(xa[khi]-xa[klo]) + (y2a[klo]/6) * \
        (-(3*a**2-1)*(xa[khi]-xa[klo])+(3*b**2-1)*y2a[khi])
    y = a*ya[klo]+b*ya[khi]+((a**3-a)*y2a[klo]+(b**3-b)*y2a[khi])*(h**2)/6
    return y, y1


# #Testing on the cubic function
# xx = range(11)
# yy = [y**3 for y in xx]
# n =  len(xx)
#
# y2 = spline(xx,yy,0,300)
#
# xs = [i*0.01 for i in range(0,1001)]
# ys = [splint(xx,yy,y2,n,x)[0] for x in xs]
#
# plt.xlim(0,10)
# plt.ylim(0,1000)
#
# plt.plot(xx,yy,marker='o',markersize=4)
#   
# plt.plot(xs,ys)
# plt.show()
if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Please provide a input file")
        exit(1)

def create_curves():
    data = load_foil(sys.argv[1])  # do a list of array with data

    X = np.linspace(0, 1, 100)
    Y2e = spline(data[1], data[2], 0, 0)
    Ye = [splint(data[1], data[2], Y2e, int(data[0][0]), x)[0] for x in X]

    Y2i = spline(data[3], data[4], 0, 0)

    Yi = [splint(data[3], data[4], Y2i, int(data[0][1]), x)[0] for x in X]

    #plt.xlim(0, 1)
    #plt.ylim(-0.3, 0.3)
    #plt.plot(X, Ye, '-g')
    #plt.plot(X, Yi, '-g')
    # plt.show()
    
    return [X, Ye, Yi]

#Cette partie retourne la courbe du profil aérodynamique de l'aile.

#[X,Ye,Yi] = create_curves()
#plt.xlim(0, 1)
#plt.ylim(-0.3, 0.3)
#plt.plot(X, Ye, '-g')
#plt.plot(X, Yi, '-g')
#plt.show()
