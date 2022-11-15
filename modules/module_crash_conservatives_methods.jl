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


function ShallowWatersWhitTopology!(F,U,Fpars,xᵢ)
    # η ≡ U[2] ; u ≡ U[1]
    x₀,x₁,p,A,g=Fpars;
    hpars=(x₀,x₁,p,A);
    h=topology_function(xᵢ,hpars)
    F[1] = -(0.5*U[1]*U[1]+g*U[2])
    F[2] = -(U[2]+h)*U[1]
end

# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCIONES AUXILIARES
function topology_function(xᵢ,par_topology)
    x₀,x₁,p,A=par_topology
    (x₀ < xᵢ < x₁) ? h=-(A*(4.0*(xᵢ-x₀)*(xᵢ-x₁)*(1.0/((x₁-x₀)^2)))^p) : h=0.0
    return h;
end

function topology_function_v2(x,par_topology)
    Z=[topology_function(x[i],par_topology) for i in 1:length(x)]
    return -Z;
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
    if (N_FIELDS == 2); @. u[:,2] = 0.0; end
    return u;
end

# shallow waters whit topology (⟺ NFIELDS==2)
function SW_oscillator_condition_NFIELDS2(N,N_FIELDS,x,CIpar)
    u0,ω,A=CIpar;
    #Inicializamos el dato
    u = Array{Float64}(undef,N,N_FIELDS);
    @. u[:,2] = u0 + A*cos(ω*x);
    @. u[:,1] = 0.0
    return u;
end

function step_function_NFIELDS1(N,N_FIELDS,x,CIpar)
    x₁,x₂,A =CIpar;
    #Inicializamos el dato
    u = Array{Float64}(undef,N,N_FIELDS);
    for i in 1:length(x)
        (x₁ <= x[i] <= x₂) ? u[i,1]=A : u[i,1]=0.0;
    end
    return u;
end

function step_function_NFIELDS1_v2(N,N_FIELDS,x,CIpar)
    x₁,x₂,A₁,A₂ =CIpar;
    #Inicializamos el dato
    u = Array{Float64}(undef,N,N_FIELDS);
    for i in 1:length(x)
        (x₁ <= x[i] <= x₂) ? u[i,1]=A₁ : u[i,1]=A₂;
    end
    return u;
end

function step_function_NFIELDS2(N,N_FIELDS,x,CIpar)
    x₁,x₂,A =CIpar;
    #Inicializamos el dato
    u = Array{Float64}(undef,N,N_FIELDS);
    for i in 1:length(x)
        (x₁ <= x[i] <= x₂) ? u[i,1]=A : u[i,1]=0.0;
    end
    @. u[:,2] = 0.0;
    return u;
end

function SincSquaredPulseFunction_NFIELDS1(N,N_FIELDS,x,CIpar)
    x₁,x₂,A=CIpar;
    a=2*π*(1.0/(x₁-x₂));
    b=(x₁+x₂)*π*(1.0/(x₁-x₂));
    #Inicializamos el dato
    u = Array{Float64}(undef,N,N_FIELDS);
    for i in 1:length(x)
        if (x₁ < x[i] < x₂)
            if (x[i] ≠ b*(1.0/a))
                u[i,1]=A*(sin(a*x[i]-b)*(1.0/(a*x[i]-b)))^2;
            else
                u[i,1]=A;
            end
        else
            u[i,1]=0.0
        end
    end
    return u;
end

# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RESOLVEMOS EDO
function resolveEDO(Problem,SpaceMethod,TimeMethod,N,x_range,tend,StateSaveat;N_FIELDS=1,θ=2.0,CIfunction=oscillator_condition,CIpar=(1.0,1.0,1.0),topologyparams)
    #=
    ---- entradas obligatorias
    Problem     → problema físico a resolver
    SpaceMethod → método numérico para integrar espacialmente
    TimeMethod  → método numérico para integrar temporalmente
    N           → cantidad de números espaciales
    x_range     → rango espacial
    tend        → tiempo final de evolución
    StateSaveat → guardar (o no) para todo Δt
    ---- entradas opcionales
    N_FIELDS    → cantidad de campos a utilizar
    θ           → parámetro específico para método KT
    CIfunction  → condicion inicial
    CIpar       → parametros de la condición inicial
    =#


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
    elseif (Problem == :ShallowWatersWhitTopology!)
        SpeedMax = shallowwatersspeed;  # nombre de función velocidad máxima
        Flux_x! = ShallowWatersWhitTopology!;        # nombre de función flujo
        eqpars = push!(topologyparams,9.81);                  # parámetro (aceleración de la gravedad)
        println("Elegida la ecuación de Shallow Waters whit Topology");
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
        (Problem == :ShallowWatersWhitTopology!) ? scheme! = mp5_ShallowWaterWhitTopology! : scheme! = mp5!
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

# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ADVECTION EQUATION: CONSTANT VLEOCITY

# Finite Difference operator whitout PERIODIC Boundary Conditions
function D2(du,u,dt)
    @. du[2:end-1] = (u[3:end] - u[1:end-2])/(2*dt) # centered 2nd order
    du[1] = (-1.5*u[1]+2.0*u[2]-0.5*u[3])/dt  # forward 2nd order
    du[end] = (1.5*u[end]-2.0*u[end-1]+0.5*u[end-2])/dt # backward 2nd order
    return du
end

function velocity_advection(sol,x,Δt,Δx)
    
    dim_row=length(sol.t)
    dim_column=length(x)

    # creamos y seteamos valores de la matriz
    u_matrix = zeros(dim_row,dim_column)
    for (t_index,t_value) in enumerate(sol.t)
        # c/fila corresponde a u(x,t) para un t fijo y ∀x
        u_matrix[t_index,:] = sol(t_value)[:]
    end

    # creamos y seteamos valores de matriz
    x_diff_matrix = zeros(dim_row,dim_column)
    for t_index in 1:dim_row
        # c/fila corresponde a ∂ₓ[u(x,t)] para un t fijo y ∀x con PBC
        x_diff_matrix[t_index,:] = D2_pbc(x_diff_matrix[t_index,:],u_matrix[t_index,:],Δx)
    end
    
    # creamos y seteamos valores de matriz
    t_diff_matrix = zeros(dim_row,dim_column)
    for x_index in 1:dim_column
        # c/fila corresponde a ∂ₜ[u(x,t)] para un x fijo y ∀t sin cond. borde
        t_diff_matrix[:,x_index] = D2(t_diff_matrix[:,x_index],u_matrix[:,x_index],Δt)
    end

    # creamos y seteamos valores de matriz
    F_matrix = zeros(dim_row,dim_column)
    # c/elemento corresponde a ∂ₜ[u(x,t)]/∂ₓ[u(x,t)]
    F_matrix .= -1.0 .* (t_diff_matrix ./ x_diff_matrix)
    return F_matrix
end

# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ANALISIS DEL CHOQUE

function FindNode(sol,xvector,Δx)

    # Searching for an element in a 1D ARRAY (array1)
    array1 = abs.(sol[:,1,1]); # vector solución a tiempo inicial
    # definimos valor para buscar dentro de array1
    ϵ = Δx; sch = 0.0 + 54*Δx;

    # creamos array sólo con elementos que cumplen con u(x,t) ≤ 54Δx ≈ 0
    elementArray = filter( x -> x <= sch, array1 )
    if (length(elementArray) == 0) # caso en que el arreglo filtrado no tenga elementos
        println("Element not found.")
    else # caso en que el arreglo filtrado tenga al menos un element
        # imprimimos el valor más cercano a cero
        u_x0=minimum(elementArray);
        println("Element found in the array.","\nu(x₀,t) = ",u_x0)
    end

    x0=0.0;x0_index=0 # POR QUE?!!!!!!
    for x_index in 1:length(xvector)
        if (abs(sol[x_index,1,1]) == u_x0)
            x0_index=x_index;x0=xvector[x0_index]
            println("x₀ = ",x0,"\nindex of x₀ = ",x0_index)
        end
    end

    return u_x0,x0,x0_index
end

# Finite Difference operator whit PBC
function D2_pbc(du,u,dx)
    @. du[2:end-1] = (u[3:end] - u[1:end-2])/(2*dx)
    du[1] = (u[2] - u[end])/(2*dx)
    du[end] = (u[1] - u[end-1])/(2*dx)
    return du
end

function x_diff_function(sol,x,Δx,Δt)
    
    dim_row=length(sol.t)
    dim_column=length(x)

    # creamos y seteamos valores de la matriz
    u_matrix = zeros(dim_row,dim_column)
    for (t_index,t_value) in enumerate(sol.t)
        # c/fila corresponde a u(x,t) para un t fijo y ∀x
        u_matrix[t_index,:] = sol(t_value)[:]
    end

    # creamos y seteamos valores de matriz
    x_diff_matrix = zeros(dim_row,dim_column)
    for t_index in 1:dim_row
        # c/fila corresponde a ∂ₓ[u(x,t)] para un t fijo y ∀x con PBC
        x_diff_matrix[t_index,:] = D2_pbc(x_diff_matrix[t_index,:],u_matrix[t_index,:],Δx)
    end

    return x_diff_matrix
end

function FindDivergence(sol,x_diff,xvector,x0_index,sch)
    # buscamos max{|∂ₓ[u(x,t)]|}
    # seteamos tupla de resultados en cero
    max_xdiff=(0.0,0.0,0.0);    # (max{|∂ₓ[u(x,t)]|},tₛ,xₛ)
    max_index=(0.0,0.0);        # (t₀Index,x₀Index)

    for (t_index,t_value) in enumerate(sol.t) # hacemos una busqueda para cada tiempo
        # esto nos devuelve una tupla con el valor y el indice del maximo encontrado
        aux_max = findmax(abs.(x_diff[t_index,:]));
        # actualizamos valor del maximo e imponemos condicion de que el nodo no se desplace
        if (aux_max[1]≥max_xdiff[1]) && (sol[x0_index,1,t_index]≤sch)
            max_xdiff = (aux_max[1],t_value,xvector[aux_max[2]])
            max_index = (t_index,aux_max[2])
        end
    end

    # calculamos el invariante = area bajo la curva en el tiempo de choque
    secction_cte=0.5*max_xdiff[3]*(sol[max_index[2],1,max_index[1]]+sol[1,1,1]);
    # calculamos la vleocidad de propagación del choque
    crash_velocity=0.5*sol[max_index[2],1,max_index[1]];

    println("(max{|∂ₓ[u(x,t)]|},tₛ,xₛ)=",max_xdiff,"\n(t₀Index,x₀Index)=",max_index);
    println("constant section = ",secction_cte,"\ncrash velocity = ",crash_velocity);
    
    return max_xdiff,max_index,secction_cte,crash_velocity
end

# ----------------------------------------------------------------