#======================================================================
Funciones que permiten resolver numéricamente ecuaciones de la forma

dy/dt = f(y, par, t)

=======================================================================#

"""
EulerStep(y, f, t, dt, par)

Método de Euler 

y valor anterior

f(y,t,par) función a integrar

t tiempo

dt paso temporal

par parámetros en f
"""
function EulerStep(y, f, t, dt, par)
    """Método de Euler"""
    return y + dt*f(y, t, par)          #Paso del método Euler
end

function RK2step(y,f,t,dt, par)
    """Método de Runge-Kutta """
    #=Inserte su código aquí=#
end

function RK4step(y,f,t,dt, par)
    """Método de Runge-Kutta 4"""
    #=Inserte su código aquí=#
end

function myODEproblem(f, y0, intervalo, par)
    """Devolvemos f, y0, el intervalo temporal, y los parámetros de la función en una tupla"""
    return (f, y0, intervalo, par)
end



function myODEsolver(Problem, Method; dt::Float64 = 0.01)
    #Evolucionamos el problem Problem utilizando el método Method con saltos demporales dt
    f, y0, intervalo, par = Problem        #Datos específicos del problema que estamos resolviendo
    tini, tfin = intervalo                 #Tiempos iniciales y finales
    N = 1 + Int(floor((tfin-tini)/dt))     #Cantidad de pasos temporales
    y = Array{typeof(y0[1])}(undef, N, length(y0))               #Vector donde guardaremos y
    t = zeros(N)                           #Vector donde guardaremos y
    y[1,:] .= y0                           #Dato inicial
    t[1] = tini
    for i in 2:N
        t[i] = tini + (i-1)*dt
        y[i,:] .= Method(y[i-1,:], f, t[i-1],dt, par)   #Pasos temporales
    end
    return (t ,y)
end
