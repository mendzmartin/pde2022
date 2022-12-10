function make_model(grid_type,params)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal",1)

    # 1D line grid
    if (grid_type == "simple_line")

        path,name,dom,MeshSize=params
        
        gmsh.model.add(name)
        
        gmsh.model.geo.addPoint(dom[1],0,0,MeshSize,1)    # 1 punto vértice izquierdo
        gmsh.model.geo.addPoint(dom[2],0,0,MeshSize,2)    # 2 punto vértice derecho

        gmsh.model.geo.addLine(1,2,1)             # linea que une puntos 1 y 2
        gmsh.model.geo.synchronize()
        
        gmsh.model.geo.addPhysicalGroup(0,[1],3) # grupo formado por punto izquierdo
        gmsh.model.setPhysicalName(0,3,"left_point")    # le damos nombre al grupo
        gmsh.model.geo.synchronize()

        gmsh.model.geo.addPhysicalGroup(0,[2],4) # grupo formado por punto derecho
        gmsh.model.setPhysicalName(0,4,"right_point")   # le damos nombre al grupo
        gmsh.model.geo.synchronize()
        
        gmsh.model.geo.addPhysicalGroup(1,[1],1) # grupo formado por linea
        gmsh.model.setPhysicalName(1,1,"segment") # le damos nombre a la linea
        gmsh.model.geo.synchronize()

        gmsh.model.mesh.generate(1)
        gmsh.write(path*name*".msh")
        gmsh.finalize()
        model=GmshDiscreteModel(path*name*".msh")
    end

    # 2D simple rectangle grid
    if (grid_type == "simple_rectangle")

        path,name,dom,MeshSize,quad_state=params

        # gmsh.option.setNumber(name, value)
        gmsh.option.setNumber("General.Terminal",0)
        gmsh.model.add(name)

        # creamos puntos vértice
        # gmsh.model.geo.addPoint(x,y,z,meshSize=0.,tag=-1)
        gmsh.model.geo.addPoint(dom[1],dom[3],0,MeshSize,1)             # 1 vertice inferior izquierdo
        gmsh.model.geo.addPoint(dom[2],dom[3],0,MeshSize,2)       # 2 vértice inferior derecho
        gmsh.model.geo.addPoint(dom[2],dom[4],0,MeshSize,3)   # 3 vértice superior derecho
        gmsh.model.geo.addPoint(dom[1],dom[4],0,MeshSize,4)        # 4 vértice superior izquierdo
        # creamos lineas de unión entre vértices
        # gmsh.model.geo.addLine(startTag,endTag,tag=-1)
        gmsh.model.geo.addLine(1,2,1) # 1 linea inferior
        gmsh.model.geo.addLine(2,3,2) # 2 línea lateral derecha
        gmsh.model.geo.addLine(3,4,3) # 3 linea superior
        gmsh.model.geo.addLine(4,1,4) # 4 linea lateral izquierda
        # creamos curva de unión entre lineas
        gmsh.model.geo.addCurveLoop([1,2,3,4],100) # the rectangle
        gmsh.model.geo.synchronize()

        # make the surface
        # gmsh.model.geo.addPlaneSurface(wireTags,tag=-1)
        gmsh.model.geo.addPlaneSurface([100],101) # the surface
        gmsh.model.geo.synchronize()

        #=
            creamos grupos para definir condiciones de bordes
            gmsh.model.geo.addPhysicalGroup(dim,tags,tag=-1,name="")
            gmsh.model.setPhysicalName(dim,tag,name)
        =#
        # creamos grupo físico de puntos vértices
        gmsh.model.geo.addPhysicalGroup(0,[1,2,3,4],300)  # grupo 0D formado por cuatro puntos
        gmsh.model.setPhysicalName(0,300,"ext_points")
        gmsh.model.geo.synchronize()

        # creamos grupo físico de curva de lineas externa
        # gmsh.model.geo.addPhysicalGroup(1,[100],301)      # de esta forma no funciona!!
        gmsh.model.geo.addPhysicalGroup(1,[1,2,3,4],301)    # grupo 1D formado por cuatro lineas
        gmsh.model.setPhysicalName(1,301,"ext_lines")
        gmsh.model.geo.synchronize()

        # creamos grupo físico con superficie interna del rectángulo
        gmsh.model.addPhysicalGroup(2,[101],302)
        gmsh.model.setPhysicalName(2,302,"surface")          # grupo 2D formado por una superficie
        gmsh.model.geo.synchronize()

        if quad_state
            # gmsh.model.mesh.setRecombine(dim,tag,angle=45.)
            gmsh.model.mesh.setRecombine(2,302) # for 2D quadrilaterals
        end

        # gmsh.model.mesh.generate(dim=3)
        gmsh.model.mesh.generate(2)
        gmsh.write(path*name*".msh")
        gmsh.finalize()
        model=GmshDiscreteModel(path*name*".msh")
    end

    return model;
end

#=
    utils links
    + https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.jl
=#