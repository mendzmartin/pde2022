#!/usr/bin/julia

"""
    RUN COMMANDS
    Via REPL => julia
                include("module_SquareDomain_NavierStokes_SteadyStateIncomp.jl")
    Via Bash => chmod +x module_SquareDomain_NavierStokes_SteadyStateIncomp.jl
                ./module_SquareDomain_NavierStokes_SteadyStateIncomp.jl
"""


using Gridap;

# import Pkg; Pkg.add("LineSearches");
using LineSearches: BackTracking;

function weak_problem_NeumannAndDirichlet(conv,dconv,dΩ,dΓ,∇ue)
    a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ;
    c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ - ∫(v⊙∇ue)dΓ;
    dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ;
    return a,c,dc;
end

function weak_problem_else(conv,dconv,dΩ)
    a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ;
    c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ;
    dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ;
    return a,c,dc;
end



function run_NavierStokes(n,dom,BCtype,Re,iter;IntFunct=lagrangian,OrderIntFunct=2,MeasureDegree=2)
#=
    n       ≝ número de elementos finitos
    dom     ≝ coordenadas que definen el dominio (x₀,x₁,y₀,y₁)
    BCtype  ≝ tipo de condiciones de contorno (eg: "original","fulldirichlet","NeumannAndDirichlet")
    Re      ≝ Número de Reynolds
    iter    ≝ número de iteraciones para integrador no lineal
    IntFunct ≝ función de interpolación (eg: lagrangian)
    OrderIntFunct ≝ orden de la función de interpolación para las velocidades
    MeasureDegree   ≝ Degree of Measure function to define dΩ
=#

    # grilla de tamaño n²
    partition = (n,n);

    # creamos el modelo con elementos cartesianos
    model = CartesianDiscreteModel(dom,partition);
    writevtk(model,"NS_steady");

    if BCtype=="original"
        println("BCtype = ", BCtype);
        labels = get_face_labeling(model);
        add_tag_from_tags!(labels,"D1",[1,2,3,4,5,7,8]); # linea inferior + lineas laterales + 4 puntos vértice
        add_tag_from_tags!(labels,"D2",[6,]);            # linea superior
    elseif BCtype=="fulldirichlet"
        println("BCtype = ", BCtype);
        # modified
        labels = get_face_labeling(model);
        add_tag_from_tags!(labels,"D1",[1,7,3,2,4,8]);   # lineas laterales + 4 puntos vértice
        add_tag_from_tags!(labels,"D2",[5,6]);           # lineas superior e inferior
    elseif BCtype=="NeumannAndDirichlet"
        println("BCtype = ", BCtype);
        ue(x)=VectorValue(1-x[2]*x[2], 0);
        ∇ₙue(x)=VectorValue(0, 0);
        labels = get_face_labeling(model);
        # dirichlet conditions
        add_tag_from_tags!(labels,"D1",[1,7,3]);         # linea lateral izquierda + 2 puntos vértice
        add_tag_from_tags!(labels,"D2",[5,6]);           # linea superior e inferior
        # neumann conditions
        add_tag_from_tags!(labels,"N1",[2,4,8]);         # linea lateral derecha + 2 puntos vértice
        Γ=BoundaryTriangulation(model,tags="N1");
        dΓ=Measure(Γ,2);
        nb=get_normal_vector(Γ);

        # https://gridap.github.io/Tutorials/dev/pages/t002_validation/
        #=
        We need to tell the Gridap library that the gradient of the function u is available in
        the function ∇u (at this moment u and ∇u are two standard Julia functions without any connection between them).
        This is done by adding an extra method to the function gradient (aka ∇) defined in Gridap
        =#
        # import Gridap: ∇;
        # ∇(::typeof(ue)) = ∇ue;
    end

    # FE SPACES
    # FE spaces for velocities
    reffe_vel = ReferenceFE(IntFunct,VectorValue{2,Float64},OrderIntFunct);
    if BCtype=="original"
        V = TestFESpace(model,reffe_vel,conformity=:H1,labels=labels,dirichlet_tags=["D1","D2"]);
    elseif ((BCtype=="fulldirichlet") || (BCtype=="NeumannAndDirichlet"))
        V = TestFESpace(model,reffe_vel,conformity=:H1,labels=labels,dirichlet_tags=["D1","D2"]);
    end

    # FE spaces for pressure
    reffe_press = ReferenceFE(IntFunct,Float64,OrderIntFunct-1;space=:P);
    Q = TestFESpace(model,reffe_press,conformity=:L2,constraint=:zeromean);

    if BCtype=="original"
        uD1 = VectorValue(0,0);
        uD2 = VectorValue(1,0);
    elseif ((BCtype=="fulldirichlet") || (BCtype=="NeumannAndDirichlet"))
        uD1(x) = VectorValue(1-x[2]*x[2],0);
        uD2 = VectorValue(0,0);
    end

    U = TrialFESpace(V,[uD1,uD2]);
    P = TrialFESpace(Q);

    Y = MultiFieldFESpace([V, Q]);
    X = MultiFieldFESpace([U, P]);

    # TRIANGULATION AND INTEGRATION QUADRATURE
    # Space and Boundary of the sistem
    Ωₕ = Triangulation(model);
    dΩ = Measure(Ωₕ,MeasureDegree);

    # NONLINEAR WEAK FORM

    conv(u,∇u) = Re*(∇u')⋅u;
    dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u);

    # re-definition from https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations
    if BCtype=="NeumannAndDirichlet"
        a,c,dc = weak_problem_NeumannAndDirichlet(conv,dconv,dΩ,dΓ,∇ₙue);
    else
        a,c,dc = weak_problem_else(conv,dconv,dΩ);
    end

    # bilineal forms
    res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v);
    jac((u,p),(du,dp),(v,q)) = a((du,dp),(v,q)) + dc(u,du,v);

    op = FEOperator(res,jac,X,Y);

    nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=iter);
    solver = FESolver(nls);

    uh, ph = solve(solver,op);

    name_folder="images/P05_IncompNavierStokesRe$(Re)_"*BCtype*"_results";
    writevtk(Ωₕ,name_folder,cellfields=["uh"=>uh,"ph"=>ph]);
end

