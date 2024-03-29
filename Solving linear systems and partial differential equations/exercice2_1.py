import numpy as np
import exercice1_2 as ex
from random import randint


"""
def gradient_conj(A, b, x):
    r = np.subtract(b, np.matmul(A, x))
    p = r
    rsold = np.matmul(r.transpose(), r)

    for i in range(1, len(b)):
        Ap = np.matmul(A, p)
        alpha = rsold/np.matmul(p.transpose(), Ap)
        x = np.append(alpha*p, x)
        r = np.subtract(alpha*Ap, r)
        rsnew = np.matmul(r.transpose(), r)
        if np.sqrt(rsnew) < 10 ^ (-10):
            break
        p = r + (rsnew/rsold) * p
        rsold = rsnew

"""


def gradient_conj(A, b, x):
    r = b - np.dot(A, x)
    p = r
    rsold = np.dot(r.transpose(), r)

    for i in range(1, len(b)):
        Ap = np.dot(A, p)
        alpha = rsold/np.dot(p.transpose(), Ap)
        print("add " + str(alpha*p) + "to x")
        x = x + alpha*p
        print("new x = " + str(x))
        r = r - alpha*Ap
        rsnew = np.dot(r.transpose(), r)
        if np.sqrt(rsnew) < 10 ^ (-10):
            break
        p = r + (rsnew/rsold) * p
        rsold = rsnew
        print("x at the end : " + str(x))
        return x


A = ex.mat_sym_creuse(5, 0)

b = (np.array([randint(1, 100) for i in range(5)])).transpose()
x = (np.array([0 for i in range(5)])).transpose()

print("x = ")
print(x)
print("A = ")
print(A)

new = gradient_conj(A, b, x)

print("new x = ")
print(new)
print(x)
