#!/usr/bin/julia

"""
    RUN COMMANDS
    Via REPL => julia
                include("module_poisson_example_01.jl")
    Via Bash => chmod +x module_poisson_example_01.jl
                ./module_poisson_example_01.jl
"""



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## paquetes necesarios

# activamos el proyecto "gridap_makie" donde se intalarán todos los paquetes
import Pkg; Pkg.activate("gridap_makie");

using Gridap;
using GridapGmsh;

# import Pkg; Pkg.add("LineSearches")
using LineSearches: BackTracking;

# Pkg.add("IterativeSolvers")
using IterativeSolvers

# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## FUNCIONES DE PRUEBA Y FUENTES

function solution_and_source_01(param)
    x₀,y₀=param;
    u(x) = 4*((x[1]-x₀)^2 - (x[2]-y₀)^3) - 5.0*x[2]; # test solution
    ∇u(x) = 8.0 - 24*(x[2] - y₀);
    f(x) = -∇u(x);                     # source
    return u,f,∇u;
end

function solution_and_source_02(param)
    x₀,y₀=param;
    u(x) = (x[1]-x₀)+(x[2]-y₀);    # test solution
    ∇u(x) = 0.0;
    f(x) = -∇u(x);                 # source
    return u,f,∇u;
end

function solution_and_source_03(param)
    x₀,y₀,a,b,p,q=param;
    u(x) = a*((x[1]-x₀)^p)+b*((x[2]-y₀)^q);                            # test solution
    ∇u(x) = a*p*(p-1)*((x[1]-x₀)^(p-2))+b*q*(q-1)*((x[2]-y₀)^(q-2));
    f(x) = -∇u(x);                                                     # source
    return u,f,∇u;
end
# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## SOLVER FUNCTIONS

#=
    Linear solver phase
        https://gridap.github.io/Tutorials/dev/pages/t001_poisson/#Solver-phase-1
=#
function LU_solver(params)
    ls = LUSolver(); solver = LinearFESolver(ls);
    return solver;
end

#=
    import Pkg; Pkg.add("LineSearches");
    using LineSearches: BackTracking;
    Non-linear sovler phase
        https://gridap.github.io/Tutorials/dev/pages/t008_inc_navier_stokes/#Nonlinear-solver-phase-1
        https://github.com/JuliaNLSolvers/NLsolve.jl
=#
function NL_solver(params)
    method,tol,iter=params
    # nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(),ftol = 1e-16, iterations=5);
    nls = NLSolver(show_trace=true, method=method, linesearch=BackTracking(),ftol = tol, iterations=iter);
    solver = FESolver(nls);
    return solver;
end
# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## DEFINITION OF ABSTACT WEAK PROBLEM

function weak_problem_full_dirichlet(dΩ,f)
    a(u,v) = ∫( ∇(u)⋅∇(v) )*dΩ; # in a(u,v) are all the dependencies with u that are unknowns
    b(v) = ∫(v*f )*dΩ;          # here everything that is source
    return a,b;
end

function weak_problem_int_dirichlet(dΩ,f,nb,ue,dΓ)
    a(u,v) = ∫( ∇(u)⋅∇(v) )*dΩ; # in a(u,v) are all the dependencies with u that are unknowns
    b(v) = ∫(v*f )*dΩ + ∫(v*(nb ⋅ ∇(ue)))*dΓ;   # here everything that is source
    return a,b;
end
# ----------------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## FUNCIÓN PRINCIPAL PARA RESOLVER POISSON CON FEM

function poisson_solver(TypeConf,MeasureDegree,BC,OrderIntFunct,TSSfunction,TSSfparams;TypeMesh="coarse",IntFunct=lagrangian,SolvFunction=LU_solver,SolvFuncPar=(0.0),CreateDir=false)

    #=
    TypeConf        ≝ Type of Configuration (eg: 00,01,etc)
    MeasureDegree   ≝ Degree of Measure function to define dΩ
    BC              ≝ Boundary Condition (eg: int_dirichlet, full_dirichlet)
    OrderIntFunct   ≝ Order of IntFunct variable
    TSSfunction     ≝ Function name to define test solution and source (eg:solution_and_source_01,etc)
    TSSfparams      ≝ Parameters to TSSfunction variable
    TypeMesh        ≝ Type of mesh (eg: "fine", "coarse")
    IntFunct        ≝ Interpolation function (eg: Lagrangian)
    SolvFunction    ≝ Name of solver function (eg:LU_solver,NL_solver)
    SolvFuncPar     ≝ Parameters to solver function
    CreateDir       ≝ Create directories (eg: false/true)
    =#


    # en caso de que no estén creados los directorios, los creamos.
    if (CreateDir==true)
        mkdir("models");
        mkdir("images");
    end

    if (TypeMesh=="fine")
        model = GmshDiscreteModel("models/rectangle_hole_fine.msh");
    elseif (TypeMesh=="coarse")
        model = GmshDiscreteModel("models/rectangle_hole_coarse.msh");
    end

    # Space and Boundary of the sistem
    Ω = Triangulation(model);
    dΩ = Measure(Ω,MeasureDegree);

    # Reference of finite elements
    reffe = ReferenceFE(IntFunct,Float64,OrderIntFunct);

    if (BC == "int_dirichlet")
        dirichlet_tags="int" ;
    elseif (BC == "full_dirichlet")
        dirichlet_tags=["int","ext"];
    end

    V0 = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags=dirichlet_tags);

    # test function and source
    ue,f=TSSfunction(TSSfparams);

    # internal Dirichlet boundary condition
    U = TrialFESpace(V0,ue);
    Γ₁ = BoundaryTriangulation(model,tags=dirichlet_tags);
    n₁ = get_normal_vector(Γ₁);

    if (BC == "int_dirichlet")
        neumann_tags = "ext"
        #neumanntags = "int"
        Γ = BoundaryTriangulation(model,tags=neumann_tags);
        dΓ = Measure(Γ,MeasureDegree);
        nb = get_normal_vector(Γ);
    end

    if (BC == "full_dirichlet")
        a,b = weak_problem_full_dirichlet(dΩ,f);
    elseif (BC == "int_dirichlet")
        a,b = weak_problem_int_dirichlet(dΩ,f,nb,ue,dΓ);
    end

    # # generate system like Ax=b with Gridap.jl
    op = AffineFEOperator(a,b,U,V0);

    # solve the problem
    uh = solve(SolvFunction(SolvFuncPar),op);

    # save solution data
    if (BC == "full_dirichlet")
        writevtk(Ω,"images/P01_conf"*TypeConf*"_solución_dir",cellfields=["uh_dir" => uh]);
    elseif (BC == "int_dirichlet")
        writevtk(Ω,"images/P01_conf"*TypeConf*"solución_newmann",cellfields=["uh_neu" => uh]);
    end

    # error between exact solution (test solution) and numerical solution
    e = ue-uh;
    # save error data
    if (BC == "full_dirichlet")
        writevtk(Ω,"images/P01_conf"*TypeConf*"_error_dir",cellfields=["e_dir" => e]);
    elseif (BC == "int_dirichlet")
        writevtk(Ω,"images/P01_conf"*TypeConf*"_error_newmann",cellfields=["e_neu" => e]);
    end

    # L²-norm and H₁-norm calculations
    el2 = sqrt(sum( ∫( e*e )*dΩ ));
    eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ));
    uh1 = sqrt(sum( ∫( uh*uh + ∇(uh)⋅∇(uh) )*dΩ ));
    println("l2 error = ",el2,"\nh1 error = ",eh1/uh1);
    
    return ue,uh,e,el2,eh1/uh1;
end
# ----------------------------------------------------------------