#!/usr/bin/julia

"""
    RUN COMMANDS
    Via REPL => julia
                include("myODEfunctions.jl")
    Via Bash => chmod +x myODEfunctions.jl
                ./myODEfunctions.jl
"""

#======================================================================
Funciones que permiten resolver numéricamente ecuaciones de la forma
    dy/dt = f(y, param, t)
=======================================================================#

""" 
    y             := valor anterior
    f(y,t,param)  := función a integrar
    t             := tiempo
    dt            := paso temporal
    param         := paramámetros en f
"""

"""Método de Euler"""
function EulerStep(y,f,t,dt,param)
    return y + dt*f(y,t,param)          #Paso del método Euler
end

"""
    Método de Runge-Kutta 2do Orden
    Método de Heun con un solo corrector (a2=1/2)
"""
function RK2step_Heun(y,f,t,dt,param)
    k_1 = f(y,t,param)                      # fuerza a tiempo no evolucionado

    y_improved = y + (k_1*dt)               # solución semi-evolucionada
    k_2 = f(y_improved,(t+dt),param)        # fuerza a tiempo evolucionado

    # y(i+1) = (y(i) + increment_function)
    y = y + 0.5*dt*(k_1 + k_2)              # solución a tiempo evolucionado
    return y
end

"""
    Método de Runge-Kutta 2do Orden
    Método del punto medio (a2=1)
"""
function RK2step_midpoint(y,f,t,dt,param)
    k_1 = f(y,t,param)                          # fuerza a tiempo no evolucionado
    dt_improved = 0.5*dt                        # medio paso temporal

    y_improved = y + (k_1*dt_improved)          # solución semi-evolucionada
    k_2 = f(y_improved,(t-dt_improved),param)   # fuerza a tiempo semievolucionado

    # y(i+1) = (y(i) + increment_function)
    y = y + dt*k_2                              # solución a tiempo evolucionado
    return y
end

"""
    Método de Runge-Kutta 2do Orden
    Método Ralston (a2 = 2/3)
"""
function RK2step_Ralston(y,f,t,dt,param)
    k_1 = f(y,t,param)                          # fuerza a tiempo no evolucionado
    dt_improved = 0.25*dt                       # medio paso temporal

    y_improved = y + (3.0*k_1*dt_improved)      # solución semi-evolucionada
    k_2 = f(y_improved,(t-dt_improved),param)   # fuerza a tiempo semievolucionado

    # y(i+1) = (y(i) + increment_function)
    y = y + (1.0/3.0)*dt*(k_1+2.0*k_2)          # solución a tiempo evolucionado
    return y
end


"""
    Método de Runge-Kutta 4to Orden
    Método Ralston (a2 = 2/3)
"""
function RK4step(y,f,t,dt, param)
    #=Inserte su código aquí=#
    k_1 = f(y,t,param)                             # fuerza a tiempo no evolucionado
    dt_improved = 0.5*dt                           # medio paso temporal

    y_improved_k2 = y + (k_1*dt_improved)          # solución semi-evolucionada
    k_2 = f(y_improved_k2,(t-dt_improved),param)   # fuerza a tiempo semievolucionado

    y_improved_k3 = y + (k_2*dt_improved)          # solución semi-evolucionada
    k_3 = f(y_improved_k3,(t-dt_improved),param)   # fuerza a tiempo semievolucionado

    y_improved_k4 = y + (k_3*dt)                   # solución semi-evolucionada
    k_4 = f(y_improved_k4,t,param)                 # fuerza a tiempo semievolucionado

    # y(i+1) = (y(i) + increment_function)
    y = y + (1.0/6.0)*dt*(k_1+2.0*(k_2+k_3)+k_4)   # solución a tiempo evolucionado
    return y
end

function myODEproblem(f, y0, intervalo, param)
    """Devolvemos f, y0, el intervalo temporal, y los paramámetros de la función en una tupla"""
    return (f, y0, intervalo, param)
end


# Evolucionamos el probleama utilizando el método con saltos temporales dt
function myODEsolver(Problem, Method; dt::Float64 = 0.01)
    f,y0,intervalo,param = Problem               # Datos específicos del problema que estamos resolviendo
    tini,tfin=intervalo                          # Tiempos iniciales y finales
    N=1+Int(floor((tfin-tini)/dt))               # Cantidad de pasos temporales
    y=Array{typeof(y0[1])}(undef,N,length(y0))   # Vector donde guardaremos y
    t=zeros(N)                                   # Vector donde guardaremos t
    y[1,:] .= y0                                 # Dato inicial
    t[1]=tini
    for i in 2:N
        t[i] = tini + (i-1)*dt
        y[i,:] .= Method(y[i-1,:], f, t[i-1],dt, param)   # Pasos temporales
    end
    return (t,y)
end