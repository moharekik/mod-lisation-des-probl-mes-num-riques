import sys
from airfoil import load_foil
from interpolate import create_curves
from interpolate import spline
from interpolate import splint
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Please provide a input file")
        exit(1)

def flux_laminaires():
    # on charge les courbes de l'aile
    X, Ye, Yi = create_curves()
    # lambda
    Lambdas = np.linspace(0, 1, 20)
    Hmin = -0.3
    Hmax = 0.3
    # plusieurs courbes des flux laminaires, extrados et intrados
    Yefls = []
    Yifls = []

    # calcule des Y
    for l in Lambdas:
        Yefl = []
        Yifl = []
        for i in range(len(Ye)):
            ye = (1 - l)*Ye[i] + l*Hmax
            yi = (1 - l)*Yi[i] + l*Hmin
            Yefl.append(ye)
            Yifl.append(yi)
        Yefls.append(Yefl)
        Yifls.append(Yifl)
        Yefl = []
        Yifl = []

    plt.xlim(0, 1)
    plt.ylim(-0.3, 0.3)
    plt.plot(X, Ye, '-g', linewidth=3)
    plt.plot(X, Yi, '-g', linewidth=3)
    for i in range(len(Lambdas)):
        plt.plot(X, Yefls[i], '-g')
        plt.plot(X, Yifls[i], '-g')
    plt.show()


            #pression = PS + (1/2)*d*Speeds[i]**2

flux_laminaires()
