#!/usr/bin/julia

"""
    RUN COMMANDS
    Via REPL => julia
                include("wave_equation_finite_differences.jl")
    Via Bash => chmod +x wave_equation_finite_differences.jl
                ./wave_equation_finite_differences.jl
"""

function create_D_2_per(N)
    D_2_per = sparse(Tridiagonal([-0.5 for i in 1:N-1],[0.0 for i in 1:N],[0.5 for i in 1:N-1]))
    D_2_per[1,end] = -0.5
    D_2_per[end,1] = 0.5
    dropzeros!(D_2_per)
    return D_2_per
end

# Operador diferencia finita de 2do orden
function create_D2_2_per(N)
    D2_2_per = BandedMatrix{Float64}(Zeros(N,N), (N-1,N-1))

    # definición de coeficientes no triviales
    a=Float64(-2);b=Float64(1);

    # seteo de diagonal principal
    D2_2_per[band(0)] .= a

    # llenado de bandas principales
    D2_2_per[band(1)] .= b
    D2_2_per[band(-1)] .= b
    
    # llenado de bandas triángulo superior
    D2_2_per[band(N-1)] .= b

    # llenado de bandas triángulo inferior
    D2_2_per[band(-N+1)] .= b
    
    D2_2_per = sparse(D2_2_per)
    dropzeros!(D2_2_per)
    return D2_2_per
end

# Operador diferencia finita de 4to orden
function create_D_4_per(N)
    D_4_per = BandedMatrix{Float64}(Zeros(N,N), (N-1,N-1))

    # definición de coeficientes no triviales
    a=Float64(-1/12);b=Float64(2/3);

    # seteo de diagonal principal
    D_4_per[band(0)] .= 0.0

    # llenado de bandas principales
    D_4_per[band(1)] .= b
    D_4_per[band(-1)] .= -b
    D_4_per[band(2)] .= a
    D_4_per[band(-2)] .= -a
    
    # llenado de bandas triángulo superior
    D_4_per[band(N-1)] .= -b
    D_4_per[band(N-2)] .= -a
    
    # llenado de bandas triángulo inferior
    D_4_per[band(-N+1)] .= b
    D_4_per[band(-N+2)] .= a
    
    D_4_per = sparse(D_4_per)
    dropzeros!(D_4_per)
    return D_4_per
end

# Operador diferencia finita de 6to orden
function create_D_6_per(N)
    D_6_per = BandedMatrix{Float64}(Zeros(N,N), (N-1,N-1))

    # definición de coeficientes no triviales
    a=Float64(1/60);b=Float64(-3/20);c=Float64(3/4);

    # seteo de diagonal principal
    D_6_per[band(0)] .= 0.0

    # llenado de bandas principales
    D_6_per[band(1)]  .= c
    D_6_per[band(-1)] .= -c
    D_6_per[band(2)]  .= b
    D_6_per[band(-2)] .= -b
    D_6_per[band(3)]  .= a
    D_6_per[band(-3)] .= -a
    
    # llenado de bandas triángulo superior
    D_6_per[band(N-1)] .= -c
    D_6_per[band(N-2)] .= -b
    D_6_per[band(N-3)] .= -a
    
    # llenado de bandas triángulo inferior
    D_6_per[band(-N+1)] .= c
    D_6_per[band(-N+2)] .= b
    D_6_per[band(-N+3)] .= a
    
    D_6_per = sparse(D_6_per)
    dropzeros!(D_6_per)
    return D_6_per
end

# Operador diferencia finita de 8vo orden
function create_D_8_per(N)
    D_8_per = BandedMatrix{Float64}(Zeros(N,N), (N-1,N-1))

    # definición de coeficientes
    a=Float64(-1/280);b=Float64(4/105);c=Float64(-1/5);d=Float64(4/5);

    # seteo de diagonal principal
    D_8_per[band(0)]    .= 0.0

    # llenado de bandas principales
    D_8_per[band(1)]    .= d
    D_8_per[band(-1)]   .= -d
    D_8_per[band(2)]    .= c
    D_8_per[band(-2)]   .= -c
    D_8_per[band(3)]    .= b
    D_8_per[band(-3)]   .= -b
    D_8_per[band(4)]    .= a
    D_8_per[band(-4)]   .= -a
    
    # llenado de bandas triángulo superior
    D_8_per[band(N-1)] .= -d
    D_8_per[band(N-2)] .= -c
    D_8_per[band(N-3)] .= -b
    D_8_per[band(N-4)] .= -a
    
    # llenado de bandas triángulo inferior
    D_8_per[band(-N+1)] .= d
    D_8_per[band(-N+2)] .= c
    D_8_per[band(-N+3)] .= b
    D_8_per[band(-N+4)] .= a
    
    D_8_per = sparse(D_8_per)
    dropzeros!(D_8_per)
    return D_8_per
end

# Definimos discretización espacial (lado derecho de la ecuación diferencial)
function F!(dr,r,p,t)
    # second order version
    dx, D = p
    h = 1. /dx
    u = @view r[:,1]
    v = @view r[:,2]
    du = @view dr[:,1]
    dv = @view dr[:,2]
    mul!(du, D, v, h, 0)  #du/dt = h*Dv
    mul!(dv, D, u, h, 0)  #dv/dt = h*Du
    #Nota: mul!(C, A, B, α, β) hace la operación α*A*B + β*C y la guarda en C
end

function Burgers!(du,u,p,t)
    # first order version
    dx, D = p
    h = 1. /dx
    #Nota: mul!(C, A, B, α, β) hace la operación α*A*B + β*C y la guarda en C
    u2 = u .^ 2
    mul!(du, D, u2, -(h*0.5), 0)  #du = -0.5*h*D(u²)
end

function Burgers_Crapodina!(du,u,p,t)
    # first order version
    dx, D = p
    h = 1. /dx
    #Nota: mul!(C, A, B, α, β) hace la operación α*A*B + β*C y la guarda en C
    #mul!(du, D, u, -h, 0)  #du = -h*D(u)
    du .= -u.*(D*u)*h
    #du=u.*du
end

function resolve_EDO_space( FinDiffOp,  # Operador de diferencia finita
                            Δx,         # Space step
                            Δt,         # Time step
                            func,       # Función del lado derecho de la EDO
                            CI,         # Condiciones iniciales
                            T,          # Tiempo final de evolución
                            Method)     # Método a usar para la resolver EDO

    time_tuple = (0.0,T)    # tupla con intervalo temporal
    param_tuple = (Δx, FinDiffOp) # tupla con parámetros
    prob = ODEProblem(func,CI,time_tuple,param_tuple);
    sol = solve(prob,Method(),dt=Δt,saveat=Δt,adaptive=false);
    return sol;
end

# cálculo de la energía para la ecuación de ondas
function energy(sol,    # solución dada por paquete DifferentialEquations
                L,      # Intervalo total espacial
                T,      # Tiempo final
                t_end)  # 

    KE=Array{Float64}(undef,(t_end+1))  # término de energía cinética
    PE=Array{Float64}(undef,(t_end+1))  # término de energía potencial
    E=Array{Float64}(undef,(t_end+1))   # energía total

    for i in 0:t_end
        KE[i+1]=sum(abs2,sol(i*T/t_end)[:,2])*0.5*L
        PE[i+1]=sum(abs2,sol(i*T/t_end)[:,1])*0.5*L
        E[i+1]=KE[i+1]+PE[i+1]
    end
    return E
end

# Factor Q
function Q(sol1, sol2, sol4, t, i)
    return norm(sol1(t)[1:end,i] - sol2(t)[1:2:end,i])/norm(sol2(t)[1:2:end,i] - sol4(t)[1:4:end,i])
end