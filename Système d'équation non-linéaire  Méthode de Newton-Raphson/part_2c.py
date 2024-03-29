import numpy as np
from part_1 import Newton_Raphson_backtrack, Newton_Raphson
import matplotlib.pyplot as plt
from numpy.polynomial import legendre


def E(X):
    s = 0
    for x in X:
        s1 = 0
        for y in X:
            if x != y:
                s1 += np.log10(abs(x-y))
        s += np.log10(abs(x+1)) + np.log10(abs(x-1)) + 1/2 * s1
    return s


def derivativeE(X, i):
    s = 0
    for xj in X:
        if xj != X[i]:
            s += 1/(X[i]-xj)
    return 1/(X[i]+1) + 1/(X[i]-1) + 1/2*s


def gradE(X):
    return [derivativeE(X, i) for i in range(len(X))]


def jacobGradE(X):
    N = len(X)
    J = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j:
                s = 0
                for xj in X:
                    if xj != X[i]:
                        s += 1/((X[i]-xj)**2)
                J[i][i] = -1/((X[i]+1)**2) - 1/((X[i]-1)**2) - 1/2*s
            else:
                J[i][j] = 1/(2*((X[i]-X[j])**2))
    return J


# Electrostatic equilibrium equation solving
N = 6
X0 = np.linspace(-0.99, 0.99, N)
S = Newton_Raphson(gradE, jacobGradE, X0, 100, 0.001)
print("Electrostatic equilibrium solution for N=6 \n", S)


# Plot
Y = np.zeros(N)
plt.plot([-1, 1], [0, 0], color='black')
plt.scatter([-1, 1], [0, 0], color='red', label='charge électrostatique')
plt.scatter(S, Y, color='blue', label='position d\'équilibre statique')

# Legendre


def legendreRoots():
    coefficients = np.zeros(N + 1)
    coefficients[-1] = 1
    leg_poly = legendre.Legendre(coefficients)
    derivative_roots = leg_poly.deriv().roots()
    gl_points = np.concatenate([[-1], derivative_roots, [1]])

    return gl_points


gl_points = legendreRoots()
print(" Roots of the derivative of the Legendre polynomials \n", gl_points)
plt.scatter(gl_points[1:N], Y[1:], color="orange",
            label="racine de la dérivé du polynome de Legendre")
plt.legend()
plt.show()
# EXTREMUM


def extremumType(S):
    eigvals = np.linalg.eigvals(jacobGradE(S))
    if all(eigvals > 0):
        print("The solution corresponds to a minimum.")
    elif all(eigvals < 0):
        print("The solution corresponds to a maximum.")
    else:
        print("The nature of the solution is unclear.")


extremumType(S)
