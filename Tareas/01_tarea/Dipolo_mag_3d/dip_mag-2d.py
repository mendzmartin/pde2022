# Computa numericamente las orbitas de un electron en un campo magnetico externo de forma
# dipolar
#
import matplotlib as mpl
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from scipy.integrate import ode
import matplotlib.animation as animation

# Vamos a integrar numericamente las siguientes ecuaciones:
# 
# (las cuales escribimos en primer orden como:)
# rho' = sigma
# sigma' = p**2/rho**3(1-mu/rho)(1-2mu/rho)
# phi' = p/rho**2(1-mu/rho)
#
# reescaleando podemos poner mu=p=1.
mu=1.
p=1.
#
def h(r,s,p,mu):
    return s*s/2. + np.power(p,2.) * np.power((1.-mu/r)/r,2.)

plt.figure()
dx = 1./400.
X = np.arange(0.01,4.01,dx)
Y = np.arange(-3.01,3.01,dx)
X, Y = np.meshgrid(X,Y)


CS = plt.contour(X, Y, h(X,Y,4.,1.),[0.,0.1,0.5,0.85,1.,1.1,1.5,3.])
#CS = plt.contour(X, Y, h(X,Y,1.,1.))


plt.clabel(CS, inline=1, fontsize=10)
plt.title('Simplest default with labels')
plt.show()
# Primero las condiciones iniciales U=(rho,sigma,phi)
N=3
x0 = 0.
rho0 = 10.
pert = 0.0000000000000001
#pert = 0.0000001
sigma0 = -p * np.sqrt(1./16./mu/mu - (1.-mu/rho0)*(1.-mu/rho0)/rho0/rho0) + pert

# la velocidad radial esta elegida para que la particula llegue al maximo del potencial y regrese
# (tardando un tiempo infinito para hacerlo) perturbandola un poquito incrementandola ingresa a la zona interna,
# con una perturbacion que disminuye la velocidad inical no alcanza a ingresar.

u0 = [rho0, sigma0, np.arctan(-p*(1.-mu/rho0)/rho0/sigma0)]

# Luego definimos la funcion del integrador, o sea la ecuacion

def f(x, u):
    return [u[1], p * p / (u[0]*u[0]*u[0]) * (1. - mu / u[0]) * (1. - 2. * mu / u[0]), p / u[0] /u[0] * (1. - mu / u[0])]

# Definimos cual integrador usar
# distintos tipos, para este problema que es muy sensible es recomendable usar la dop853, pero pruebe tambien la dopri5

#r = ode(f).set_integrator('zvode', method='bdf')
#r = ode(f).set_integrator('dopri5')  # Runge Kutta de orden 4/5
r = ode(f).set_integrator('dop853')   # Runge Kutta de orden 8
r.rtol = 1.e-16

# Damos los valores iniciales

r.set_initial_value(u0, x0)

# El valor final de integracion

x1 = 10000 # con estos datos iniciales llegamos a y=0 para este valor de x

# Esto es solo a los fines de graficar, si da un error: out of bounds poner K mas grande
dx = 0.1
K=100000
#dx=1.
#K=10000
KK=0
U = np.zeros((K,N))
X = np.zeros(K)

# Integramos K veces separadas (para graficar)
i = 0
while r.successful() and r.t < x1 and r.y[0] < rho0*1.2:
    r.integrate(r.t+dx)
    U[i,:] = r.y
    X[i] = r.t
    i=i+1
    KK=KK+1

# Finalmente graficamos la solucion    

ax = plt.subplot(111, projection='polar')
#ax.plot(U[:,2],U[:,0], '.', color='r', linewidth=1)
ax.plot(U[:,2],U[:,0], color='r', linewidth=1)
ax.set_rmax(6.10)
ax.grid(True)

ax.set_title("z=0 orbits", va='bottom')
plt.show()

#plt.plot(X[:],U[:,0])
#plt.show()
    
# MATE EL GRAFICO PARA CONTINUAR!


fig = plt.figure()
ax1 = fig.add_subplot(111, projection='polar' )
#ax1.set_rmax(6.10)
ax1.set_rmax(10.)
ax1.grid()
line, = ax1.plot([], [], lw=1)
time_template = 'time = %.1fs'
time_text1 = ax1.text(0.01, 0.95, '', transform=ax1.transAxes)


def init():
    line.set_data([], [])
    ax1.set_rmax(10.)
    time_text1.set_text('')
    return line, time_text1

def animate(i):
    line.set_data(U[0:i,2],U[0:i,0])
    time_text1.set_text(time_template % (i*dx))
  #   time_text2.set_text(time_template % (i*dt))
    return line, time_text1

ani = animation.FuncAnimation(fig, animate, KK,
                              interval=1, blit=True, init_func=init, repeat=False)

#ani.save('double_pendulum.mp4', fps=15)
plt.show()
exit()
