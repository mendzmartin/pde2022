#!/usr/bin/julia

"""
    RUN COMMANDS
    Via REPL => julia
                include("module_crash_conservatives_methods.jl")
    Via Bash => chmod +x module_crash_conservatives_methods.jl
                ./module_crash_conservatives_methods.jl
"""

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## VELOCIDADES MÁXIMAS DE PROPAGACIÓN

function advectionspeed(U, c)
    return abs(c);
end

function burgersspeed(U, c)
    return maximum(abs, U);  #no encuentro forma de escribir esto sin que aloque memoria...
end

function shallowwatersspeed(U,c)
    return maximum(abs,U)
end

# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FLUJOS

function advection!(F,U,c)
    @. F = c*U;
end

function burgers!(F,U,Fpars)
    @. F = 0.5*U*U;
end

function shallowwaters(F,U,g)
    F[1] = -U[2]
    F[2] = -(U[2]*U[2]/U[1] + g*U[1]*U[1]/2)
end

# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## CONDICIONES INICIALES

function oscillator_condition(N,N_FIELDS,x,CIpar)
    u0,ω,A=CIpar;
    #Inicializamos el dato
    u = Array{Float64}(undef,N,N_FIELDS);
    @. u[:,1] = u0 + A*sin(ω*x);
    return u;
end

function SW_oscillator_condition(N,N_FIELDS,x,CIpar)
    u0,ω,A=CIpar;
    #Inicializamos el dato
    u = Array{Float64}(undef,N,N_FIELDS);
    @. u[:,1] = u0 + A*sin(ω*x);
    @. u[:,2] = 0.0;
    return u;
end

# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RESOLVEMOS EDO
function resolveEDO(Problem,SpaceMethod,TimeMethod,N,x_range,tend,StateSaveat;N_FIELDS=1,θ=2.0,CIfunction=oscillator_condition,CIpar=(1.0,1.0,1.0))
    start=x_range[1];stop=x_range[2];
    x=range(start,stop=stop,length=N+1)[1:end-1]; # de manera que no incluya el último punto
    Δx=Float64(x.step);
    h=1.0/Δx;

    u=CIfunction(N,N_FIELDS,x,CIpar);

    #Definimos el intervalo de integración y el paso Δt
    tspan = (0.0, tend);

    #CFL = Δt/Δx
    CFL=0.1;Δt=Δx*CFL;

    # Definimos parámetros según el problema físico a resolver
    if (Problem == :advection)
        SpeedMax = advectionspeed;   # nombre de función velocidad máxima
        Flux_x! = advection!;        # nombre de función flujo
        eqpars = 1.0;                # parámetro (velocidad de advección)
        println("Elegida la ecuación de Advección");
    elseif (Problem == :burgers)
        SpeedMax = burgersspeed;     # nombre de función velocidad máxima
        Flux_x! = burgers!;          # nombre de función flujo
        eqpars = false;              # parámetro
        println("Elegida la ecuación de Burgers");
    elseif (Problem == :shallowwaters)
        SpeedMax = shallowwatersspeed;  # nombre de función velocidad máxima
        Flux_x! = shallowwaters;        # nombre de función flujo
        eqpars = 9.81;                  # parámetro (aceleración de la gravedad)
        println("Elegida la ecuación de Shallow Waters");
    end

    # Definimos parámetros según el método  de resolución espacial elegido
    if (SpaceMethod == :KurganovTadmor) # Métod de Kurganov Tadmor
        # θ ≝ Este valor tiene que estar entre 1 y 2. Mientras más cerca de 2, menor disipación.
        # Para sistemas de ecuaciones es mejor que esté más cerca de 1 para evitar oscilaciones.

        auxvectors = createKTauxvectors(N_FIELDS);
        scheme! = KT!;
        par = (eqpars,h,θ,Flux_x!,SpeedMax,N,N_FIELDS,auxvectors);
        println("Elegido el método KT");
    elseif (SpaceMethod == :MonotonicityPreserving5) # Método de Monotonicity Preserving 5
        auxvectors = createMP5auxvectors(N_FIELDS);
        scheme! = mp5!;
        par = (eqpars,h,N,N_FIELDS,Flux_x!,SpeedMax,auxvectors);
        println("Elegido el método MP5");
    end

    # definimos el problema
    prob = ODEProblem(scheme!,u,tspan,par);
    # resolvemos con método específico para integración temporal
    if (All_saveat == true)
        sol = solve(prob,TimeMethod(),dt=Δt,saveat=Δt); # Esto es un método TVD
    else (All_saveat == false)
        sol = solve(prob,TimeMethod(),dt=Δt,saveat=t_end/100); # Esto es un método TVD
    end

    return sol,x,Δx,Δt;
end