#!/usr/bin/julia

#=
    RUN COMMANDS
    Via REPL => julia
                include("module_schrodinger_equation_testing.jl")
    Via Bash => chmod +x module_schrodinger_equation_testing.jl
                ./module_schrodinger_equation_testing.jl
=#

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Definimos rutas a directorios específicos para buscar o guardar datos
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

path_models         = "../models/";
path_images         = "../images/";
path_modules        = "../modules/"
path_gridap_makie   = "../gridap_makie/";
path_videos         = "./videos/";
path_plots          = "./plots/";


#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Activamos proyecto e intalamos paquetes para FEM
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

# activamos el proyecto "gridap_makie" donde se intalarán todos los paquetes
import Pkg; Pkg.activate(path_gridap_makie);

install_packages=false;
if install_packages
    import Pkg
    Pkg.add("GridapGmsh");
    Pkg.add("Gmsh");
    Pkg.add("FileIO");
end

using Gridap;
using GridapGmsh;
using Gmsh;
using Gridap.CellData;  # para construir condición inicial interpolando una función conocida
using Gridap.FESpaces;  # para crear matrices afines a partir de formas bilineales
using Gridap.Algebra;   # para utilizar operaciones algebraicas con Gridap
# using Gridap.Arrays
# using Gridap.ReferenceFEs

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

include(path_modules*"module_eigen_prototype.jl");  # módulo para resolver problema de autovalores
include(path_models*"mesh_generator.jl"); # módulo para construir grilla (1D)

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Seteo de variables globales
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

# declaramos parámetros constantes
const m=1.0;                            # masa
const ω=1.0;                            # frecuencia angular
const ħ=1.0;                            # constante de Planck
const x₁=0.0;                           # posición donde se centra el 1er osc. armónico
const x₂=2.0;                           # posición donde se centra el 2do osc. armónico
const γ=0.1;                            # constante de acoplamiento
const α=im*ħ*0.5*(1.0/m);               # factor multiplicativo energía cinética
const αconst=-im*0.5*m*(ω*ω)*(1.0/ħ);   # factor multiplicativo potencial armónico

@printf("VARIABLES GLOBALES:\n");
@printf("m=%.4f (mass)\nω=%.4f (frecuency)\nħ=%.4f (Planck constant)\nγ=%.4f (coupling)\n",m,ω,ħ,γ);
@printf("x₁=%.4f x₂=%.4f (QHO origin position)\n",x₁,x₂);

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
qₕ_2D(x) = 0.5*m*(ω*ω)*((x[1]-x₁)*(x[1]-x₁)+(x[2]-y₁)*(x[2]-y₁));   # oscilador armónico 2D centrado en (x₁,y₁)
rₕ(x) = 1.0;

# Formas bilineales para problema de autovalores (espacios complejos)
#  deben verificar la integración por partes
function bilineal_forms(p,q,r,dΩ)
    a(u,v) = ∫(p*∇(v)⋅∇(u)+q*v*u)*dΩ;
    b(u,v) = ∫(r*u*v)dΩ;
    return a,b;
end

# Formas bilineales para problema de autovalores (espacios real e imaginario separados)
# function bilineal_forms_ReImParts(p,q,r,dΩ)
#     a((u₁,u₂),(v₁,v₂)) = ∫(p*(∇(v₁)⋅∇(u₁)+∇(v₂)⋅∇(u₂))+q*(v₁*u₁+v₂*u₂))dΩ;
#     b((u₁,u₂),(v₁,v₂)) = ∫(r*(v₁*u₁+v₂*u₂))dΩ;
#     return a,b;
# end

function bilineal_forms_ReImParts(p,q,r,dΩ)
    a₁((u₁,v₁))=∫(p*(∇(v₁)⋅∇(u₁))+q*(v₁*u₁))dΩ;
    b₁((u₁,v₁))=∫(r*(v₁*u₁))dΩ;

    a₂((u₂,v₂))=∫(p*(∇(v₂)⋅∇(u₂))+q*(v₂*u₂))dΩ;
    b₂((u₂,v₂))=∫(r*(v₂*u₂))dΩ;

    a((u₁,u₂),(v₁,v₂)) = a₁((u₁,v₁))+a₂((u₂,v₂))
    b((u₁,u₂),(v₁,v₂)) = b₁((u₁,v₁))+b₂((u₂,v₂))
    return a,b;
end


# Norma L₂
function norm_L2(u,dΩ)
    return sqrt(real(sum(∫(u'*u)*dΩ)));
end

# funciones para hamiltoniano 2x2 1D
α₁(x)=αconst*(x[1]-x₁)*(x[1]-x₁); # oscilador armónico 1D centrado en x₁
α₂(x)=αconst*(x[1]-x₂)*(x[1]-x₂); # oscilador armónico 1D centrado en x₂

#= +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++ Funciones útiles para el problema de autovalores completo
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

# funciones para problema de autovalores (Ec. de Sturm Liouville)
pH(x) = 0.5*(ħ*ħ)*(1.0/m);                                          # factor para energía cinética
qH₁(x) = 0.5*m*(ω*ω)*(x[1]-x₁)*(x[1]-x₁);                           # oscilador armónico 1D centrado en x₁
qH₂(x) = 0.5*m*(ω*ω)*(x[1]-x₂)*(x[1]-x₂);                           # oscilador armónico 1D centrado en x₂
rH(x) = 1.0;
sH(x) = γ;

function bilineal_forms_eigenprob_H(p,q₁,q₂,r,s,dΩ)
    a((u₁,u₂),(v₁,v₂)) = ∫(p*(∇(v₁)⋅∇(u₁)+∇(v₂)⋅∇(u₂))+q₁*v₁*u₁+q₂*v₂*u₂+s*(v₁*u₁+v₂*u₂))*dΩ;
    b((u₁,u₂),(v₁,v₂)) = ∫(r*(v₁*u₁+v₂*u₂))dΩ;
    return a,b;
end

# function bilineal_forms_eigenprob_H_ReImParts(p,q₁,q₂,r,s,dΩ)

#     # parte real de la 1er coordenada
#     a₁((u₁,u₃),(v₁,v₃))=∫(p*(∇(v₁)⋅∇(u₁))+q₁*(v₁*u₁)+s*(v₃*u₃))*dΩ;
#     b₁((u₁,v₁))=∫(r*(v₁*u₁))*dΩ;

#     # parte imaginaria de la 1er coordenada
#     a₂((u₂,u₄),(v₂,v₄))=∫(p*(∇(v₂)⋅∇(u₂))+q₁*(v₂*u₂)+s*(v₄*u₄))*dΩ;
#     b₂((u₂,v₂))=∫(r*(v₂*u₂))*dΩ;

#     # parte real de la 2da coordenada
#     a₃((u₃,u₁),(v₃,v₁))=∫(p*(∇(v₃)⋅∇(u₃))+q₂*(v₃*u₃)+s*(v₁*u₁))*dΩ;
#     b₃((u₃,v₃))=∫(r*(v₃*u₃))*dΩ;

#     # parte imaginaria de la 2da coordenada
#     a₄((u₄,u₂),(v₄,v₂))=∫(p*(∇(v₄)⋅∇(u₄))+q₂*(v₄*u₄)+s*(v₂*u₂))*dΩ;
#     b₄((u₄,v₄))=∫(r*(v₄*u₄))*dΩ;

#     a((u₁,u₂,u₃,u₄),(v₁,v₂,v₃,v₄)) = a₁((u₁,u₃),(v₁,v₃))+a₂((u₂,u₄),(v₂,v₄))+a₃((u₃,u₁),(v₃,v₁))+a₄((u₄,u₂),(v₄,v₂))
#     b((u₁,u₂,u₃,u₄),(v₁,v₂,v₃,v₄)) = b₁((u₁,v₁))+b₂((u₂,v₂))+b₃((u₃,v₃))+b₄((u₄,v₄))

#     # a((u₁,u₂,u₃,u₄),(v₁,v₂,v₃,v₄)) = ∫(p*(∇(v₁)⋅∇(u₁)+∇(v₂)⋅∇(u₂)+∇(v₃)⋅∇(u₃)+∇(v₄)⋅∇(u₄))+q₁*(v₁*u₁+v₂*u₂)+q₂*(v₃*u₃+v₄*u₄)+s*(v₁*u₁+v₂*u₂+v₃*u₃+v₄*u₄))*dΩ;
#     # b((u₁,u₂,u₃,u₄),(v₁,v₂,v₃,v₄)) = ∫(r*(v₁*u₁+v₂*u₂+v₃*u₃+v₄*u₄))dΩ;
#     return a,b;
# end

function bilineal_forms_eigenprob_H_ReImParts(p,q₁,q₂,r,s,dΩ)

    # parte real de la 1er coordenada
    a₁((u₁,u₃),v₁)=∫(p*(∇(v₁)⋅∇(u₁))+q₁*(v₁*u₁)+s*(v₁*u₃))*dΩ;
    b₁((u₁,v₁))=∫(r*(v₁*u₁))*dΩ;

    # parte imaginaria de la 1er coordenada
    a₂((u₂,u₄),v₂)=∫(p*(∇(v₂)⋅∇(u₂))+q₁*(v₂*u₂)+s*(v₂*u₄))*dΩ;
    b₂((u₂,v₂))=∫(r*(v₂*u₂))*dΩ;

    # parte real de la 2da coordenada
    a₃((u₃,u₁),v₃)=∫(p*(∇(v₃)⋅∇(u₃))+q₂*(v₃*u₃)+s*(v₃*u₁))*dΩ;
    b₃((u₃,v₃))=∫(r*(v₃*u₃))*dΩ;

    # parte imaginaria de la 2da coordenada
    a₄((u₄,u₂),v₄)=∫(p*(∇(v₄)⋅∇(u₄))+q₂*(v₄*u₄)+s*(v₄*u₂))*dΩ;
    b₄((u₄,v₄))=∫(r*(v₄*u₄))*dΩ;

    a((u₁,u₂,u₃,u₄),(v₁,v₂,v₃,v₄)) = a₁((u₁,u₃),v₁)+a₂((u₂,u₄),v₂)+a₃((u₃,u₁),v₃)+a₄((u₄,u₂),v₄)
    b((u₁,u₂,u₃,u₄),(v₁,v₂,v₃,v₄)) = b₁((u₁,v₁))+b₂((u₂,v₂))+b₃((u₃,v₃))+b₄((u₄,v₄))

    return a,b;
end