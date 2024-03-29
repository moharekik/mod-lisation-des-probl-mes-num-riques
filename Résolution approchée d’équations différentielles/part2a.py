from part1 import meth_epsilon, step_rk4, step_euler,meth_n_step
import matplotlib.pyplot as plt
import numpy as np


gamma = 0.1 #taux de croissance intrinsèque de croissance des lions par anné
kappa = 1000 # capacité maximal du mileu (résérve de faune)


"""

Titre de l'étude :
Analyse de la dynamique de la population de lions dans une réserve de faune africaine

Introduction :
Cette étude vise à analyser la dynamique de la population de lions dans une
réserve de faune en Afrique. Nous utiliserons un modèle de population
prédateur-proie pour examiner les interactions entre les lions et leurs proies,
ainsi que l'impact de divers facteurs environnementaux sur la taille de la
population de lions.

"""

def mathus(t,y):
	return gamma*y

def verhulst(t,y):
	return gamma*y*(1-y/kappa)


def approx_period(data,step):
    # Recherche de l'indice du premier maximum
    max_index = np.argmax(data)
    
    # Recherche de l'indice du premier minimum après le premier maximum
    min_index = max_index + np.argmin(data[max_index:])
    
    # Calcul de la période approximative
    period = min_index * 2* step
    
    return period

#Résolution des equations différentielle

eps = 1e-2
t0 = 0
tf = 100

y0 = 100 #population initiale

S1mathus = meth_epsilon(t0,y0,tf,eps,mathus,step_rk4)
S1verhulst = meth_epsilon(t0,y0,tf,eps,verhulst,step_rk4)

plt.plot(np.linspace(t0,tf,len(S1mathus)),S1mathus, label ='Mathus exemple 1')
plt.plot(np.linspace(t0,tf,len(S1verhulst)),S1verhulst, label = " Verhulst exemple 1 ")

y0 = 50 #population initiale

S2mathus = meth_epsilon(t0,y0,tf,eps,mathus,step_rk4)
S2verhulst = meth_epsilon(t0,y0,tf,eps,verhulst,step_rk4)

plt.plot(np.linspace(t0,tf,len(S2mathus)),S2mathus, label ='Mathus exemple 2')
plt.plot(np.linspace(t0,tf,len(S2verhulst)),S2verhulst, label = " Verhulst exemple 2 ")

y0 = 400 #population initiale

S3mathus = meth_epsilon(t0,y0,tf,eps,mathus,step_rk4)
S3verhulst = meth_epsilon(t0,y0,tf,eps,verhulst,step_rk4)

#plt.plot(np.linspace(t0,tf,len(S3mathus)),S3mathus)
#plt.plot(np.linspace(t0,tf,len(S3verhulst)),S3verhulst)
plt.legend()
plt.show()

a = 2/3
b = 4/3
c = 0.4
d = 0.2
def lotka_volterra(t,Y):
	dNdt = Y[0] * (a - b*Y[1])
	dPdt = Y[1] * (c * Y[0] - d)
	return np.array([dNdt, dPdt])


Y0 = np.array([1.05,1.05])

S1= meth_epsilon(t0, Y0, tf, eps, lotka_volterra,step_euler)
plt.plot(np.linspace(t0,tf,len(S1)),S1)



Y0 = np.array([1.1,1.1])

S2 = meth_epsilon(t0, Y0, tf, eps, lotka_volterra,step_euler)
#plt.plot(np.linspace(t0,tf,len(S2)),S2)
plt.legend()
plt.show()

C = S1.T
periode = approx_period(C[1],(tf-t0)/len(C[1]))
print(periode)
plt.plot(C[0],C[1])
plt.show()
