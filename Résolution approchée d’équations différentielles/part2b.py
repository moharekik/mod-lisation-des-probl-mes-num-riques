import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq

##Modelisation d'un pendule à un unique maillon
# parametres
g = 9.81
L = 1.
omega0 = np.sqrt(g/L)
# fonctions de base
def F(Y,t):
	"""second membre de l'EDO"""
	global omega0
	dY = np.array([Y[1],-omega0**2*np.sin(Y[0])])
	return dY


# integration avec une CI theta0
def solution(theta0):
	T = np.linspace(0.,2*(2*np.pi)/omega0,100)
	# condition initiale
	Y0 = np.array([theta0,0.])
	# integration
	Y = odeint(F, Y0, T)
	return T,Y

# tracer de la solution 
def trace(T,Y):
	Theta=Y[:,0]
	
	# trace theta(t)
	plt.plot(T,Theta,label="$\\theta$")
	plt.title('Angle $\\theta(t)$')
	plt.xlabel('t')
	plt.legend()
	plt.show()
	return  
# calcul de la période de la solution numérique 
# determination des 2 intervalles k1-1,k1  et k2-1,k2
# correspondant au passage par zero de theta
def periode(T,Theta):
	k1=0
	while Theta[k1] > 0 : 
		k1=k1+1
	k2=k1
	while Theta[k2]<0: 
		k2=k2+1
	# interpolation lineaire
	t1 = T[k1-1] - (T[k1]-T[k1-1])*Theta[k1-1]/(Theta[k1]-Theta[k1-1])
	t2 = T[k2-1] - (T[k2]-T[k2-1])*Theta[k2-1]/(Theta[k2]-Theta[k2-1])
	# d'ou la periode
	return 2*(t2-t1)

# tracer la frequence en fonction de theta0
def tracer_freq(List):
	F=[]
	for theta in List:    
		T,Y=solution(theta)
		F.append(1/periode(T,Y[:,0]))        
	plt.plot(List,F,label="Frequence")
	f0=np.sqrt(g/L)/(2*np.pi)
	print(f0)
	plt.plot(List,[f0]*len(List),label="f0")
	plt.xlabel('"$\\theta$0"')
	plt.legend()
	plt.show()
	return 

tracer_freq(np.linspace(0.,np.pi/2,100))


##modelisation d'un pendule à deux maillons 
# Paramètres
L1 = 1.0  # Longueur du premier maillon
L2 = 0.8  # Longueur du deuxième maillon
m1 = 1.0  # Masse du premier maillon
m2 = 0.5  # Masse du deuxième maillon

# Fonctions de base
def F(Y, t):
	"""Second membre de l'EDO"""
	theta1, omega1, theta2, omega2 = Y

	# Calcul des dérivées
	dtheta1 = omega1
	domega1 = (-g * (2 * m1 + m2) * np.sin(theta1) - m2 * g * np.sin(theta1 - 2 * theta2) - 2 * np.sin(theta1 - theta2) * m2 * (omega2**2 * L2 + omega1**2 * L1 * np.cos(theta1 - theta2))) / (L1 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))

	dtheta2 = omega2
	domega2 = (2 * np.sin(theta1 - theta2) * (omega1**2 * L1 * (m1 + m2) + g * (m1 + m2) * np.cos(theta1) + omega2**2 * L2 * m2 * np.cos(theta1 - theta2))) / (L2 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))

	return [dtheta1, domega1, dtheta2, domega2]


# Integration avec les conditions initiales theta1_0, omega1_0, theta2_0, omega2_0
def solution(theta1_0, omega1_0, theta2_0, omega2_0):
	T = np.linspace(0., 10., 1000)  # Intervalle de temps
	Y0 = [theta1_0, omega1_0, theta2_0, omega2_0]  # Conditions initiales
	Y = odeint(F, Y0, T)  # Intégration des équations différentielles

	return T, Y


# Tracer de la solution
def trace(T, Y):
	theta1 = Y[:, 0]
	theta2 = Y[:, 2]

	# Tracer theta1(t) et theta2(t)
	plt.plot(T, theta1, label=r"$\theta_1$")
	plt.plot(T, theta2, label=r"$\theta_2$")
	plt.title('Angles $\\theta_1(t)$ et $\\theta_2(t)$')
	plt.xlabel('t')
	plt.legend()
	plt.show()


# Exemple d'utilisation
T, Y = solution(np.pi / 4, 0, np.pi / 2, 0)  # Conditions initiales : theta1 = pi/4, theta2 = pi/2
trace(T, Y)

# Exemple d'utilisation
T, Y = solution(np.pi / 4, 0, np.pi / 2, 0)  # Conditions initiales : theta1 = pi/4, theta2 = pi/2

# Calcul des coordonnées de l'extrémité du pendule
def pendule_extremite(theta1, theta2):
	x = L1 * np.sin(theta1) + L2 * np.sin(theta2)
	y = -L1 * np.cos(theta1) - L2 * np.cos(theta2)
	return x, y

x_ext, y_ext = pendule_extremite(Y[:, 0], Y[:, 2])

# Tracer de la trajectoire de l'extrémité du pendule
plt.plot(x_ext, y_ext)
plt.title("Trajectoire de l'extrémité du pendule à deux maillions")
plt.show()

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Function to calculate the derivatives of the system
def derivs(t, y, m1, m2, l1, l2, g):
    th1, th2, p1, p2 = y
    
    c, s = np.cos(th1-th2), np.sin(th1-th2)
    denom = (m1+m2)*l1 - m2*l1*c**2
    
    th1dot = p1 / (m1*l1**2) + p2 / (m1*l1**2) * c / denom
    th2dot = p2 / (m2*l2**2) - p1 / (m2*l2**2) * c / denom
    
    p1dot = -m2*l1*l2*th2dot**2*s - (m1+m2)*g*l1*np.sin(th1)
    p2dot = m2*l1*l2*th1dot**2*s - m2*g*l2*np.sin(th2)
    
    return th1dot, th2dot, p1dot, p2dot


# Initial conditions
th1_0, th2_0 = np.pi/2, np.pi/2
p1_0, p2_0 = 0, 0
y0 = [th1_0, th2_0, p1_0, p2_0]

# Time points for integration
t_start, t_end, dt = 0, 20, 0.01
t_span = np.arange(t_start, t_end, dt)

# Solve the system using solve_ivp
sol = solve_ivp(derivs, [t_start, t_end], y0, args=(m1, m2, L1, L2, g), t_eval=t_span)

# Extract the solutions for the angles
th1, th2 = sol.y[0], sol.y[1]

# Calculate the x and y coordinates of the pendulum ends
x1, y1 = L1*np.sin(th1), -L1*np.cos(th1)
x2, y2 = x1 + L2*np.sin(th2), y1 - L2*np.cos(th2)

# Plot the trajectory of the pendulum
fig, ax = plt.subplots(figsize=(5,5))
ax.set_aspect('equal')
ax.plot(x2, y2, color='blue', alpha=0.7)
ax.plot(x2[0], y2[0], marker='o', color='black')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()
