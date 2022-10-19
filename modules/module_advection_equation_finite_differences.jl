#!/usr/bin/julia

"""
    RUN COMMANDS
    Via REPL => julia
                include("module_advection_equation_finite_differences.jl")
    Via Bash => chmod +x module_advection_equation_finite_differences.jl
                ./module_advection_equation_finite_differences.jl
"""

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCIÓN PARA CREAD OPERADOR DIFERENCIA FINITA DE 2DO ORDEN
#     EN CONJUNTO CON MAPA h2 PARA CUMPLIR SBP
function create_D_2_sbp(N)
    D_2_sbp = Tridiagonal([-0.5 for i in 1:N-1],[0.0 for i in 1:N],[0.5 for i in 1:N-1])
    D_2_sbp[1,1] = -1.0
    D_2_sbp[1,2] = 1.0
    D_2_sbp[end,end] = 1.0
    D_2_sbp[end,end-1] = -1.0
    h2 = Diagonal([1. for i in 1:N])
    h2[1,1] = 0.5
    h2[N,N] = 0.5
    return h2, D_2_sbp
end

# OPERADORES DE 4to ORDEN

# Operador diferencia finita de 4to orden
function create_D4SBP(N)
    D4SBP = BandedMatrix{Float64}(Zeros(N,N), (N-1,N-1))

    # definición de coeficientes no triviales
    a=Float64(-1/12);b=Float64(2/3);

    # seteo de diagonal principal
    D4SBP[band(0)] .= 0.0

    # llenado de bandas principales
    D4SBP[band(1)] .= b
    D4SBP[band(-1)] .= -b
    D4SBP[band(2)] .= a
    D4SBP[band(-2)] .= -a
    
    # CONDICIONES DE CONTORNO
    # llenado de subespacio cercano a x=0 y a x=L
    Head=[  -1.4117647059   1.7352941176    -0.23529411765  -0.088235294118     0.0             0.0
            -0.5            0.0             0.5             0.0                 0.0             0.0
            0.093023255814  -0.68604651163  0.0             0.68604651163       -0.093023255814 0.0
            0.030612244898  0.0             -0.60204081633  0.0                 0.65306122449   -0.081632653061]

    D4SBP[1:4,1:6] .= Head
    for i=1:4,j=1:6;D4SBP[N-(i-1),N-(j-1)]=-Head[i,j];end
    
    # definición de matriz tipo sparse
    D4SBP = sparse(D4SBP)
    dropzeros!(D4SBP)
    return D4SBP
end

# mapa de 4to orden para satisfacer SBP
function h4!(N) 
    h4 = Diagonal([1. for i in 1:N])
    h4[1,1] = 17/48
    h4[2,2] = 59/48
    h4[3,3] = 43/48
    h4[4,4] = 49/48
    h4[end,end] = h4[1,1]
    h4[end-1,end-1] = h4[2,2]
    h4[end-2,end-2] = h4[3,3]
    h4[end-3,end-3] = h4[4,4]
    return h4
end
# -------------------------------------------------------------------

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCIONES QUE EVOLUCIONAN LA FUNCIÓN DE ADVECCIÓN

# para una onda que viaja hacia la derecha
function F_sbp_matricial_right!(du,u,p,t)
    # second order version
    h11, D, c,g, parg,dx = p
    h = 1. /dx
    mul!(du, D, u, -c*h,0)
    du[1] += c * (h/h11)*(g(t, parg) - u[1])
end

# para una onda que viaja hacia la izquierda
function F_sbp_matricial_left!(du,u,p,t)
    # second order version
    hNN, D, c,g, parg,dx = p
    h = 1. /dx
    mul!(du, D, u, -c*h,0)
    du[N] += -c * (h/hNN)*(g(t, parg) - u[N])
end
# -------------------------------------------------------------------

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCIONES PARA CREAR CONDICIONES DE CONTORNO
function step_function(t,ω)
    return sign(sin(t*ω))
end

function sine_squared_pulse_function(t,par)
    t < pi ? sin(t) : 0.0
end

function wave_modulation_function(t,ω)
    return sin(t*ω)*cos(2*ω*t)
end

function reflexion(t,χ) # reflection
    return -χ
end

function refraction_left(t,param)
    c1,c2,u_minus_L,u_plus_R = param
    return (2*c2*u_plus_R+(c1-c2)*u_minus_L)/(c1+c2)
end

function refraction_right(t,param)
    c1,c2,u_minus_L,u_plus_R = param
    return (2*c1*u_minus_L+(c2-c1)*u_plus_R)/(c1+c2)
end

# -------------------------------------------------------------------

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCIÓN PARA RESOLVER EDO DE ADVECCIÓN
# resolve_EDO(FinDiffOpSBP,mapSBP,Δx,c,func01,func02,param_func02,CI,tend,Method)
function resolve_EDO(FinDiffOpSBP,  # Operador de diferencia finita con SBP
    mapSBP,              # elemento de la forma/mapa para cumplir SBP
    Δx,                  # Space step
    c,                   # wave velocity
    func01,              # función para evolucionar ecuación de advección
    func02,param_func02, # función y parametros para definir condicion de borde
    CI,                  # Condiciones iniciales
    tend,                # Tiempo final de evolución
    Method)              # Método a usar para la resolver EDO

    time_tuple = (0.0,tend);    # tupla con intervalo temporal
    Δt = Δx/abs(c);             # time step
    param_tuple = (mapSBP,FinDiffOpSBP,c,func02,param_func02,Δx); # parameters
    prob = ODEProblem(func01,CI,time_tuple,param_tuple);    # problem definition
    @time sol = solve(prob,Method(),dt=Δt,adaptive=false);  # solution
    return sol;
end

function resolve_EDO_oneStep(FinDiffOpSBP,  # Operador de diferencia finita con SBP
    mapSBP,              # elemento de la forma/mapa para cumplir SBP
    Δx,                  # Space step
    c,                   # wave velocity
    func01,              # función para evolucionar ecuación de advección
    func02,param_func02, # función y parametros para definir condicion de borde
    CI,                  # Condiciones iniciales
    Method)              # Método a usar para la resolver EDO

    time_tuple = (0,1.0);    # tupla con intervalo temporal
    Δt = Δx/abs(c);             # time step
    param_tuple = (mapSBP,FinDiffOpSBP,c,func02,param_func02,Δx); # parameters
    prob = ODEProblem(func01,CI,time_tuple,param_tuple);    # problem definition
    integrator = init(prob,Method(),dt=Δt,adaptive=false)
    step!(integrator)
    return integrator;
end

function resolve_EDO_oneStep2(FinDiffOpSBP,  # Operador de diferencia finita con SBP
    mapSBP,              # elemento de la forma/mapa para cumplir SBP
    Δx,                  # Space step
    tevol,
    Δtspecif,
    c,                   # wave velocity
    func01,              # función para evolucionar ecuación de advección
    func02,param_func02, # función y parametros para definir condicion de borde
    CI,                  # Condiciones iniciales
    Method)              # Método a usar para la resolver EDO

    #Δt = Δx/abs(c);             # time step
    Δt = Δtspecif
    time_tuple = (tevol,tevol+2*Δt);    # tupla con intervalo temporal
    param_tuple = (mapSBP,FinDiffOpSBP,c,func02,param_func02,Δx); # parameters
    prob = ODEProblem(func01,CI,time_tuple,param_tuple);    # problem definition
    integrator = init(prob,Method(),dt=Δt,adaptive=false)
    step!(integrator)
    return integrator;
end
# -------------------------------------------------------------------
