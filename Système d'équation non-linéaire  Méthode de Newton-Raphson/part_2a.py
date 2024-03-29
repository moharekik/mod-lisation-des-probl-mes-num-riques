import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle 
import part_1 as prt1

def f_elast(X0, k):
    def force(X):
        return np.array([-k * (X[0] - X0[0]), -k * (X[1] - X0[1])])

    def jacob(X):
        return np.array([[-k, 0], [0, -k]])

    return force, jacob



def f_centri(X0, k):
    def force(X):
        return np.array([k * (X[0] - X0[0]), k * (X[1] - X0[1])])

    def jacob(X):
        return np.array([[k, 0], [0, k]])

    return force, jacob



def f_gravi(X0, k):
    def force(X):
        if X[0] == X0[0] and X[1] == X0[1]:
            return 0
        return np.array([-k*(X[0]-X0[0])/((X[0]-X0[0])**2+(X[1]-X0[1])**2)**(3/2), -k*(X[1]-X0[1])/(((X[0]-X0[0])**2+(X[1]-X0[1])**2)**(3/2))])

    def jacob(X):
        if X[0] == X0[0] and X[1] == X0[1]:
            return 0
        r = k/(((X[0]-X0[0])**2 + (X[1]-X0[1])**2)**(5/2))
        return r * np.array([[2*(X[0]-X0[0])**2-(X[1]-X0[1])**2, 3*(X[0]-X0[0])*(X[1]-X0[1])], [3*(X[0]-X0[0])*(X[1]-X0[1]), 2*(X[1]-X0[1])**2-(X[0]-X0[0])**2]])
    
    return force, jacob


def f(X):
    return gravi_1[0](X) + gravi_2[0](X) + centri[0](X)

def J(X):
    return gravi_1[1](X) + gravi_2[1](X) + centri[1](X)


if __name__ == "__main__":
    gravi_1 = f_gravi(np.array([0, 0]), 1)
    gravi_2 = f_gravi(np.array([1, 0]), 0.01)
    centri = f_centri(np.array([(0.01/1.01), 0]), 1) # 0.01/1.01 selon x est le barycentre des deux masses

    print("Force totale au point [1.5, 0] :")
    print(f(np.array([1.5, 0])))
    f_expct = np.array([1.00565457, 0]) # valeur donnée sur la page du cours
    print(f_expct, "attendue\n-------------")

    print("Matrice jacobienne appliquée au point [1.5, 0] :")
    print(J(np.array([1.5, 0])))
    J_expct = np.array([[1.75259259, 0], [0, 0.6237037]]) # valeur donnée sur la page du cours
    print(J_expct, "attendue\n-------------")

    Xs = []
    for i in np.arange(-10, 10, 1):
        for j in np.arange(-10, 10, 1):
            sol = prt1.Newton_Raphson(f, J, np.array([i, j]), 100, 0.01)
            if np.linalg.norm(sol) > 0.01:
                Xs.append(sol)

    Xs_unique = []

    print("Points d'équilibres du système :")
    for x in Xs:
        if not any(np.isclose(x, y, atol=0.001).all() for y in Xs_unique):
            Xs_unique.append(x)
            print(x)
    print("-------------")

    xs = [x[0] for x in Xs_unique]
    ys = [x[1] for x in Xs_unique]

    plt.scatter(xs, ys)

    for i, point in enumerate(Xs_unique):
        plt.text(point[0] + 0.1, point[1], f"P{i+1}", fontsize=12, ha="center", va="center")
    
    plt.axis('equal')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
