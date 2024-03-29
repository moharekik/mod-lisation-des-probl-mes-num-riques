import matplotlib.pyplot as plt
import numpy as np 

def step_Heun(y,t,h,f):
	k1 = f(t,y)
	k2 = f(t + h/2, y + h/2 * k1)
	return y + h*k2


def step_euler(y, t, h, f):
	return y + h*f(t,y)

def step_midpoint(y, t, h, f):
	k1 = f(t, y)
	k2 = f(t + h / 2, y + k1 * h / 2)
	return y + h * k2

def step_rk4(y,t,h,f):
	k1 = f(t,y)
	k2 = f(t +h/2,y + (h*k1)/2)
	k3 = f(t +h/2 , y +(h*k2)/2)
	k4 = f(t + h,y +h*k3)
	return y + (h/6)*(k1 +2*k2 + 2*k3 + k4)





# def meth_n_step(y0, t0, N, h, f, meth):
# 	"""
# 	Calcule la solution d'un problème de Cauchy avec une méthode donnée
# 	à partir de la condition initiale y0 à l'instant t0,
# 	sur N pas de temps de taille h.
# 	Retourne un tableau de taille N+1 avec les valeurs de y aux différents
# 	pas de temps.
# 	"""
# 	y = [y0]  # Initialise le tableau de valeurs avec y0
# 	t = t0
# 	for i in range(N):
# 		y_n = y[-1]
# 		t += h
# 		y.append(meth(y_n, t, h, f))  # Ajoute la nouvelle valeur à la fin du tableau
# 	return y

def meth_n_step(y0, t0, N, h, f, meth):
    t=t0
    y=y0
    pts_t=[t0]
    pts_y=[y0]
    for k in range (1,N+1):
        y = meth(y, t, h, f)
        pts_t.append(t)
        pts_y.append(y)
    return np.array(pts_y)

# def meth_epsilon(t0,y0,tf,eps,f,meth):
# 	Nmax = 100
# 	N = 1
# 	h = (tf-t0)
# 	Yn = np.array(meth_n_step(y0,t0,N,h,f,meth))
# 	N = 2
# 	h = h/2
# 	Y2n = np.array(meth_n_step(y0,t0,N,h,f,meth))
# 	while (N <Nmax):
# 		if np.linalg.norm(Yn - np.array(Y2n[::2]),ord =np.inf) < eps :
# 			break
# 		Yn = np.copy(Y2n)
# 		N *= 2
# 		h /= 2
# 		Y2n = np.array(meth_n_step(y0,t0,N,h,f,meth))
# 	return Y2n

def norm_epsilon(y1 ,y2, epsilon):
    if (isinstance(y1[0],float)):
        return (np.linalg.norm(y1-y2, np.inf) >=epsilon)
    else :
        n = len(y1)
        norm = []
        for i in range(n):
            if (np.linalg.norm(y1[i] - y2[i], np.inf) >= epsilon):
                return True
        return False
    
def meth_epsilon(t0, y0, tf, eps, f, meth):
    Nmax = 400
    N = 1
    h = (tf-t0)/N
    y2n = meth_n_step(y0, t0, 2*N, h/2, f, meth)
    yn = meth_n_step(y0, t0, N, h, f, meth)
    y2nnorm = np.array([y2n[i*2] for i in range(0, int((len(y2n)+1)/2))])
    while norm_epsilon(yn, y2nnorm, eps) and N < Nmax:
        N *= 2
        h = (tf-t0)/N
        yn = y2n
        y2n = meth_n_step(y0, t0, 2*N, h/2, f, meth)
        y2nnorm = np.array([y2n[i*2] for i in range(0, int((len(y2n)+1)/2))])
    return yn



def plot_tangent_field(f, xlim, ylim, n_points=20):
	x = np.linspace(xlim[0], xlim[1], n_points)
	y = np.linspace(ylim[0], ylim[1], n_points)
	X, Y = np.meshgrid(x, y)

	U = np.zeros_like(X)
	V = np.zeros_like(Y)

	for i in range(n_points):
		for j in range(n_points):
			vec = f(0, np.array([X[i, j], Y[i, j]]))
			U[i, j] = vec[0]
			V[i, j] = vec[1]

	plt.quiver(X, Y, U, V)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.show()


	
def f1(t, y):
	return y/(1+t**2)

def f2(t, y):
	return np.array([-y[1], y[0]])
	
if __name__ == "__main__":

	y0_1 = 1
	t0_1 = 0
	tf_1 = 5
	eps = 1e-5

	sol1_rk4 = meth_epsilon(t0_1, y0_1, tf_1, eps, f1, step_rk4)
	sol1_Euler = meth_epsilon(t0_1, y0_1, tf_1, eps, f1, step_euler)
	sol1_Heun = meth_epsilon(t0_1, y0_1, tf_1, eps, f1, step_Heun)
	sol1_Midpoint = meth_epsilon(t0_1, y0_1, tf_1, eps, f1, step_midpoint)


	# Pour la première équation différentielle
	t1 = np.linspace(t0_1, tf_1, len(sol1_rk4))
	exact_sol1 = np.exp(np.arctan(t1))
	plt.plot(t1, exact_sol1, label='Solution exacte')
	plt.plot(t1, sol1_rk4, label='Solution numérique - RK4')
	plt.plot(t1, sol1_Heun, label='Solution numérique - Heun')
	plt.plot(t1, sol1_Midpoint, label='Solution numérique - Midpoint')
	plt.plot(t1, sol1_Euler, label='Solution numérique - Euler')
	plt.plot(t1, exact_sol1, label='Solution exacte')
	plt.xlabel('t')
	plt.ylabel('y')
	plt.legend()
	plt.title("Comparaison des solutions pour l'équation différentielle de dimension 1")
	plt.show()


	y0_2 = np.array([1, 0])
	t0_2 = 0
	tf_2 = 5

	sol2_euler =meth_epsilon(t0_2, y0_2, tf_2, eps, f2, step_euler)
	sol2_Heun = meth_epsilon(t0_2, y0_2, tf_2, eps, f2, step_Heun)
	sol2_Midpoint = meth_epsilon(t0_2, y0_2, tf_2, eps, f2, step_midpoint)
	t2 = np.linspace(t0_2, tf_2, len(sol2_euler))
	exact_sol2 = np.array([np.cos(t2), np.sin(t2)]).T

	plt.plot(t2, [y[0] for y in sol2_Heun], label='Solution numérique y1(t) - Heun')
	plt.plot(t2, [y[1] for y in sol2_Heun], label='Solution numérique y2(t) - Heun')
	plt.plot(t2, exact_sol2[:, 0], label='Solution exacte y1(t)')
	plt.plot(t2, exact_sol2[:, 1], label='Solution exacte y2(t)')
	plt.plot(t2, [y[0] for y in sol2_Midpoint], label='Solution numérique y1(t) - Midpoint')
	plt.plot(t2, [y[1] for y in sol2_Midpoint], label='Solution numérique y2(t) - Midpoint')
	plt.plot(t2, [y[0] for y in sol2_euler], label='Solution numérique y1(t) - Euler')
	plt.plot(t2, [y[1] for y in sol2_euler], label='Solution numérique y2(t) - Euler')
	plt.plot(t2, exact_sol2[:, 0], label='Solution exacte y1(t)')
	plt.plot(t2, exact_sol2[:, 1], label='Solution exacte y2(t)')
	plt.xlabel('t')
	plt.ylabel('y')
	plt.title("Comparaison des solutions pour l'équation différentielle de dimension 2")
	plt.show()

	plot_tangent_field(f2, [-2, 2], [-2, 2])
