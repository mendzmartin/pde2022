# Computa numericamente las orbitas de un electron en un campo magnetico externo de forma
# dipolar
#
import matplotlib as mpl
import numpy as np
import matplotlib
#matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from scipy.integrate import ode
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Vamos a integrar numericamente las siguientes ecuaciones:
# 
# (las cuales escribimos en primer orden como:)
# rho' = sigma
# sigma' = p**2/rho**3(1-mu/rho)(1-2mu/rho)
# phi' = p/rho**2(1-mu/rho)
#
# reescaleando podemos poner mu=p=1.
p=-0.1 # momento angular (reescaleado con gamma y m)
v2=0.5 # cuadrado de la velocidad total
gamma = np.sqrt(1.-v2) # factor de Lorentz
a = 0.1 # amplitud del campo (reescaleado con e, gamma, m y c)
#
# Primero las condiciones iniciales U=(rho,rho_t,theta, theta_t,phi)
N=5
x0 = 0.
rho0 = 1.
rho_t0 = -0.02
#pert = 0.000000000000000
pert = 0.01
theta0 = np.pi/2.0
theta_t0 = pert
phi0 = 0

h0 = (p/rho0/np.sin(theta0) + a*np.sin(theta0)/rho0/rho0)
h20 = h0*h0
v2 = (rho_t0*rho_t0 + rho0*rho0*(theta_t0*theta_t0) + h20) / (1.+h20)

if (v2 > 0.99):
                                print('velocidad demasiado alta! v2=%f', v2)
                                exit()
gamma = np.sqrt(1.-v2)                        

u0 = [rho0, rho_t0, theta0, theta_t0, phi0]

# Luego definimos la funcion del integrador, o sea la ecuacion

def f(x, u):
    sint = np.sin(u[2])
    cost = np.cos(u[2])
    phi_t = (p/u[0]/sint + a/u[0]/u[0]*sint)/u[0]/sint/gamma
    
    return [u[1], a/u[0]/u[0]*sint*sint*phi_t/gamma + (v2 - u[1]*u[1])/u[0],u[3],
            (-2.*a/u[0]/u[0]/u[0]/gamma + phi_t)*phi_t*sint*cost - 2.*u[1]*u[3]/u[0],phi_t]

# Definimos cual integrador usar
# distintos tipos, para este problema que es muy sensible es recomendable usar la dop853, pero pruebe tambien la dopri5

#r = ode(f).set_integrator('zvode', method='bdf')
#r = ode(f).set_integrator('dopri5')  # Runge Kutta de orden 4/5
r = ode(f).set_integrator('dop853')   # Runge Kutta de orden 8
r.rtol = 1.e-16

# Damos los valores iniciales

r.set_initial_value(u0, x0)

# El valor final de integracion

#x1 = 10000 # con estos datos iniciales llegamos a y=0 para este valor de x
x1 = 5000.
# Esto es solo a los fines de graficar, si da un error: out of bounds poner K mas grande
dx = 0.1
K=50001
#K=100000
#dx=1.
#K=10000
KK=0
U = np.zeros((K,N))
X = np.zeros(K)

# Integramos K veces separadas (para graficar)
i = 0
while r.successful() and r.t < x1 and r.y[0] < rho0*12.:
    r.integrate(r.t+dx)
    U[i,:] = r.y
    X[i] = r.t
    i=i+1
    KK=KK+1

# Finalmente graficamos la solucion    

#ax = plt.subplot(111, projection='polar')
ax = plt.subplot(111, projection='3d')
#ax.plot(U[:,2],U[:,0], '.', color='r', linewidth=1)
#ax.plot(U[:,4],U[:,0], color='r', linewidth=1)
ax.plot(U[0:KK,0]*np.sin(U[0:KK,4])*np.sin(U[0:KK,2]),
        U[0:KK,0]*np.cos(U[0:KK,4])*np.sin(U[0:KK,2]),
        U[0:KK,0]*np.cos(U[0:KK,2]), color='r', linewidth=1)
#ax.set_rmax(6.10)
ax.grid(True)

ax.set_title("generic orbits", va='bottom')
plt.show()

#plt.plot(X[:],U[:,0])
#plt.show()
    
# MATE EL GRAFICO PARA CONTINUAR!


fig = plt.figure()
#ax1 = fig.add_subplot(111, projection='polar' )
#ax1.set_rmax(6.10)
#ax1.set_rmax(10.)
#ax1 = fig.add_subplot(111, projection='3d')
ax1 = fig.gca(projection='3d')
ax1.grid(True)
line, = ax1.plot([], [], [], color='r' ,lw=1)
time_template = 'time = %.1fs'
time_text1 = ax1.text(0.01, 0.01, 0.95, '', transform=ax1.transAxes)
ax1.set_xlim((-1.5, 1.5))
ax1.set_ylim((-1.5, 1.5))
ax1.set_zlim((-1.5, 1.5))

def init():
    line.set_data([], [])
    line.set_3d_properties([])
#    ax1.set_rmax(10.)
    time_text1.set_text('')
    return line, time_text1

def animate(i):
    line.set_data(U[0:i,0]*np.sin(U[0:i,4])*np.sin(U[0:i,2]),
        U[0:i,0]*np.cos(U[0:i,4])*np.sin(U[0:i,2]))
    line.set_3d_properties(U[0:i,0]*np.cos(U[0:i,2]))
    time_text1.set_text(time_template % (i*dx))
  #   time_text2.set_text(time_template % (i*dt))
    return line, time_text1

ani = animation.FuncAnimation(fig, animate, KK,
                              interval=1, blit=True, init_func=init, repeat=False)

#ani.save('dipole-3d.mp4', fps=15)
plt.show()
exit()
