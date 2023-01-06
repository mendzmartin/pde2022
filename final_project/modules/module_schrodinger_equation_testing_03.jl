#!/usr/bin/julia

#=
    RUN COMMANDS
    Via REPL => julia
                include("module_schrodinger_equation_testing_03.jl")
    Via Bash => chmod +x module_schrodinger_equation_testing_03.jl
                ./module_schrodinger_equation_testing_03.jl
=#

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Definimos rutas a directorios específicos para buscar o guardar datos
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

path_models         = "../outputs/Output_Testing_03_SingleEigenProblemAndImplicitMethod/models/";
path_images         = "../outputs/Output_Testing_03_SingleEigenProblemAndImplicitMethod/images/";
path_modules        = "../modules/"
path_gridap_makie   = "../gridap_makie/";
path_videos         = "./videos/";
path_plots          = "../outputs/Output_Testing_03_SingleEigenProblemAndImplicitMethod/plots/";

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Activamos proyecto e intalamos paquetes para FEM
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

# activamos el proyecto "gridap_makie" donde se intalarán todos los paquetes
import Pkg; Pkg.activate(path_gridap_makie);

# necesario cambiar a true si se corre por 1ra vez
install_packages=false;
if install_packages
    import Pkg
    Pkg.add("Gridap");
    Pkg.add("GridapGmsh");
    Pkg.add("Gmsh");
    Pkg.add("FileIO");
end

using Gridap;
using GridapGmsh;
using Gmsh;
using Gridap.CellData; # para construir condición inicial interpolando una función conocida
using Gridap.FESpaces; # para crear matrices afines a partir de formas bilineales
using Gridap.Algebra

using Plots;

# crear directorios en caso de no haberlo hecho
create_directories = false;
if (create_directories==true)
    mkdir(path_models);
    mkdir(path_images);
end

using FileIO;

# en caso de querer plotear dentro de Jupiter Notebook
#  debemos usar algunos paquetes. (no funciona en VSCode)
plot_s = false;
if plot_s
    using GridapMakie, GLMakie; # Para graficar 
    using FileIO;               # Gráficos y salidas
end

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Instalamos otros paquetes útiles
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#
using Printf; # para imprimir salidas con formatos

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Instalamos paquetes para operaciones algebraicas
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#
# necesario cambiar a true si se corre por 1ra vez
install_packages=false;
if install_packages
    import Pkg
    Pkg.add("LinearAlgebra");
    Pkg.add("SparseArrays");
    Pkg.add("LinearAlgebra");
    Pkg.add("Arpack");
end
using LinearAlgebra;
using SparseArrays;
using SuiteSparse;
using Arpack;

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Importamos módulos
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

include(path_modules*"module_eigen.jl");  # módulo para resolver problema de autovalores
include(path_modules*"module_mesh_generator.jl"); # módulo para construir grilla (1D)

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Seteo de variables globales
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#
# declaramos parámetros constantes
const m=1.0;const ω=1.0;const ħ=1.0;const x₁=0.0;const x₂=2.0;const γ=0.1;
const α=im*ħ*0.5*(1.0/m);const αconst=-im*0.5*m*(ω*ω)*(1.0/ħ);const β=-im*γ*(1.0/ħ);

@printf("VARIABLES GLOBALES:\n");
@printf("m=%.4f (mass)\nω=%.4f (frecuency)\nħ=%.4f (Planck constant)\nγ=%.4f (coupling)\n",m,ω,ħ,γ);
@printf("x₁=%.4f x₂=%.4ff (QHO origin position)\n",x₁,x₂);

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Funciones útiles
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

# the triangulation and integration aproximated Lebesgue measure
function measures(model,degree,tags_boundary)
    # triangulation of the integration domain
    Ω=Triangulation(model);
    dΩ=Measure(Ω,degree);
    # triangulation of the boundary domain whit boundary conditions
    Γ=BoundaryTriangulation(model,tags=tags_boundary);
    dΓ=Measure(Γ,degree)
    return Ω,dΩ,Γ,dΓ;
end
# definimos espacios de referencia
function reference_FEspaces(method,type,order)
    reff=ReferenceFE(method,type,order);
    return reff;
end

# funciones para problema de autovalores (Ec. de Sturm Liouville)
pₕ(x) = 0.5*(ħ*ħ)*(1.0/m);                                          # factor para energía cinética
qₕ(x) = 0.5*m*(ω*ω)*(x[1]-x₁)*(x[1]-x₁);                            # oscilador armónico 1D centrado en x₁
rₕ(x) = 1.0;

# Formas bilineales para problema de autovalores
#  deben verificar la integración por partes
function bilineal_forms(p,q,r,dΩ)
    a(u,v) = ∫(p*(∇(v)⋅∇(u))+q*(v*u))dΩ;
    b(u,v) = ∫(r*(u*v))dΩ;
    return a,b;
end

# Norma L₂
function norm_L2(u,dΩ)
    return sqrt(real(sum(∫(u'*u)*dΩ)));
end

# funciones para hamiltoniano 2x2 1D
α₁(x)=αconst*(x[1]-x₁)*(x[1]-x₁); # oscilador armónico 1D centrado en x₁
α₂(x)=αconst*(x[1]-x₂)*(x[1]-x₂); # oscilador armónico 1D centrado en x₂

# Formas bilineales para problema de autovalores

#  deben verificar la integración por partes
function a_bilineal_forms_2D(α₁,α₂,Δt,dΩ)
    a₁((u₁,u₂),v₁)=∫(2*(u₁*v₁)-(-α*(∇(v₁)⋅∇(u₁))+α₁*(u₁*v₁)+β*(u₂*v₁))*Δt)dΩ
    a₂((u₂,u₁),v₂)=∫(2*(u₂*v₂)-(-α*(∇(v₂)⋅∇(u₂))+α₂*(u₂*v₂)+β*(u₁*v₂))*Δt)dΩ
    a((u₁,u₂),(v₁,v₂))=a₁((u₁,u₂),v₁)+a₂((u₂,u₁),v₂)
    return a;
end

function b_bilineal_form_2D(α₁,α₂,u₀₁,u₀₂,Δt,dΩ)
    b₁(v₁)=∫(2*(u₀₁*v₁)+(-α*(∇(v₁)⋅∇(u₀₁))+α₁*(u₀₁*v₁)+β*(u₀₂*v₁))*Δt)dΩ
    b₂(v₂)=∫(2*(u₀₂*v₂)+(-α*(∇(v₂)⋅∇(u₀₂))+α₂*(u₀₂*v₂)+β*(u₀₁*v₂))*Δt)dΩ
    b((v₁,v₂))=b₁(v₁)+b₂(v₂)
    return b;
end