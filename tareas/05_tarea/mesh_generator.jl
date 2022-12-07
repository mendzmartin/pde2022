function make_model(grid_type,p)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)


    if grid_type == "test" || grid_type == "rectangle"  # simple square
        
        name, side_x, side_y, lc = p
        #first we build the rectangular boundary: 
        #gmsh.initialize()
        #gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("$name")
        #lc = 1e-2
        #lc = 1e-1
        lc_x = lc
        lc_y = lc*side_y/side_x
        gmsh.model.geo.addPoint(0, 0, 0, lc_x, 1)
        gmsh.model.geo.addPoint(side_x, 0,  0, lc_x, 2)
        gmsh.model.geo.addPoint(side_x, side_y, 0, lc_y, 3)
        gmsh.model.geo.addPoint(0, side_y, 0, lc_y, 4)

        # make the square boundary
        gmsh.model.geo.addLine(1, 2, 1)
        gmsh.model.geo.addLine(2, 3, 2)
        gmsh.model.geo.addLine(3, 4, 3)
        gmsh.model.geo.addLine(4, 1, 4)

        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 10) #the rectangle
        gmsh.model.geo.synchronize()

        #gmsh.model.geo.addPhysicalGroup(1, [10], 11 )
        gmsh.model.geo.addPhysicalGroup(1, [1, 2, 3, 4], 11 )
        gmsh.model.setPhysicalName(1, 11, "ext")
        gmsh.model.geo.synchronize()

        # make the surface

        gmsh.model.geo.addPlaneSurface([10], 100) #the surface
        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(2, [100], 101)
        gmsh.model.setPhysicalName(2, 101, "surface")
        gmsh.model.geo.synchronize()

        gmsh.model.mesh.generate(2)
        gmsh.write("models/$name.msh")
        gmsh.finalize()
        model = GmshDiscreteModel("models/$name.msh")
 
    end

    #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Bloque de código agregado para crear un nuevo test que permita
        utilizar elementos cuadriláteros y mallas estructuradas.
       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#
    if grid_type == "test_crapodina"   # simple square Quad and triangle elements

        println("Choose test crapodina 😃");
        
        # parámetros de entrada
        name, side_x, side_y, lc, numNodesHE , quad_state,structured_mesh,bumpFactor = p
        # Asignamos número de nodos en los bordes horizontal y vertical
        numNodesHE_hor,numNodesHE_ver = numNodesHE

        if (quad_state == true)
            gmsh.option.setNumber("Mesh.Algorithm", 5) # delquad
            gmsh.option.setNumber("Mesh.RecombineAll", 1)
        end
        gmsh.model.add("$name")

        # first we build the rectangular boundary
        lc_x = lc
        lc_y = lc*side_y/side_x
        gmsh.model.occ.addPoint(0, 0, 0, lc_x, 1)               # 1 vertice inferior izq
        gmsh.model.occ.addPoint(side_x, 0,  0, lc_x, 2)         # 2 vértice inferior der
        gmsh.model.occ.addPoint(side_x, side_y, 0, lc_y, 3)     # 3 vértice superior der
        gmsh.model.occ.addPoint(0, side_y, 0, lc_y, 4)          # 4 vértice superior izq

        # make the square boundary
        gmsh.model.occ.addLine(1, 2, 1) # 1 linea inferior
        gmsh.model.occ.addLine(2, 3, 2) # 2 linea lateral der 
        gmsh.model.occ.addLine(3, 4, 3) # 3 linea superior
        gmsh.model.occ.addLine(4, 1, 4) # 4 linea lateral izq

        gmsh.model.occ.addCurveLoop([1, 2, 3, 4], 10) #the rectangle
        gmsh.model.occ.synchronize()

        # make the surface
        gmsh.model.occ.addPlaneSurface([10], 100) #the surface
        gmsh.model.occ.synchronize()
        
        type_structured_mesh="AlternateLeft";

        if (structured_mesh == true)
            #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++
            The `setTransfiniteCurve()' meshing constraints explicitly specifies the
            location of the nodes on the curve.
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++ =# 
            # Creamos curvas de interpolación inferior y superior
            # bumpFactor = 0.20 (default)
            gmsh.model.mesh.setTransfiniteCurve(1, numNodesHE_hor, "Bump", bumpFactor)
            gmsh.model.mesh.setTransfiniteCurve(3, numNodesHE_hor, "Bump", bumpFactor)
            # Creamos curvas de interpolación lateral izquierda y derecha
            gmsh.model.mesh.setTransfiniteCurve(2, numNodesHE_ver, "Bump", bumpFactor)
            gmsh.model.mesh.setTransfiniteCurve(4, numNodesHE_ver, "Bump", bumpFactor)
            
            #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++
            The `setTransfiniteSurface()' meshing constraint uses a transfinite
            interpolation algorithm in the parametric plane of the surface to connect
            the nodes on the boundary using a structured grid. If the surface has more
            than 4 corner points, the corners of the transfinite interpolation have to
            be specified by hand:
            The way triangles are generated can be controlled by specifying "Left",
            "Right" or "Alternate" in `setTransfiniteSurface()' command.
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#
            # creamos malla estructurada en la cara 2D
            if (type_structured_mesh == "Default")
                gmsh.model.mesh.setTransfiniteSurface(100) # for structured mesh
            elseif (type_structured_mesh == "Left" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"Left")
            elseif (type_structured_mesh == "Right" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"Right")
            elseif (type_structured_mesh == "Alternate" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"Alternate")
            elseif (type_structured_mesh == "AlternateLeft" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"AlternateLeft")
            elseif (type_structured_mesh == "AlternateRight" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"AlternateRight")
            end
        end

        if (quad_state == true)
            println("Choose FE-quadrilaterals");
            gmsh.model.mesh.setRecombine(2, 100) # for 2D quadrilaterals
        else
            println("Choose FE-triangles (default)");
        end

        # creamos grupos para definir condiciones de bordes
        # gmsh.model.addPhysicalGroup(dimensión,elementos,tag)

        gmsh.model.geo.addPhysicalGroup(0, [1, 2, 3, 4], 12 )   # grupo formado por las 4 lineas del borde
        gmsh.model.setPhysicalName(0, 12, "ext_vertices")       # le damos nombre al grupo (dim = 0)
        gmsh.model.geo.synchronize()                            # sincronizamos para que sea visible

        gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4], 11 )       # grupo formado por las 4 lineas del borde
        gmsh.model.setPhysicalName(1, 11, "ext")                # le damos nombre al grupo (dim = 1)
        gmsh.model.occ.synchronize()                            # sincronizamos para que sea visible

        gmsh.model.addPhysicalGroup(2, [100], 101)              # grupo formado por cara 2D
        gmsh.model.setPhysicalName(2, 101, "surface")           # sincronizamos para que sea visible
        gmsh.model.occ.synchronize()                            # sincronizamos para que sea visible

        # generamos mesh 2D
        gmsh.model.mesh.generate(2)
        gmsh.write("models/$name.msh")

        # # esto abre una consola interactiva con gmsh
        # if !("-nopopup" in ARGS)
        #     gmsh.fltk.run()
        # end

        # finalizamos armado de grilla
        gmsh.finalize()

        # guardamos mesh en una variable
        model = GmshDiscreteModel("models/$name.msh")
    end

    #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#


    if grid_type == "square_circle" # square - circle
        
        name, side_x, side_y, cy_center_x, cy_center_y, cy_radious, lc = p

        gmsh.model.add(name)
        #lc = 1e-2
        #lc = 1e-1
        gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
        gmsh.model.geo.addPoint(side_x, 0,  0, lc, 2)
        gmsh.model.geo.addPoint(side_x, side_y, 0, lc, 3)
        gmsh.model.geo.addPoint(0, side_y, 0, lc, 4)

        # make the square boundary
        gmsh.model.geo.addLine(1, 2, 1)
        gmsh.model.geo.addLine(2, 3, 2)
        gmsh.model.geo.addLine(3, 4, 3)
        gmsh.model.geo.addLine(4, 1, 4)

        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 10) #the rectangle
        gmsh.model.geo.synchronize()

        #gmsh.model.geo.addPhysicalGroup(1, [10], 11 )
        gmsh.model.geo.addPhysicalGroup(1, [1, 2, 3, 4], 11 )
        gmsh.model.setPhysicalName(1, 11, "ext")
        gmsh.model.geo.synchronize()

        # add the circle

        lc_f = lc/4;
        gmsh.model.geo.addPoint(cy_center_x, cy_center_y, 0, lc_f, 5)
        gmsh.model.geo.addPoint(cy_center_x + cy_radious, cy_center_y, 0, lc_f, 6)
        gmsh.model.geo.addPoint(cy_center_x , cy_center_y + cy_radious, 0, lc_f, 7)
        gmsh.model.geo.addPoint(cy_center_x - cy_radious, cy_center_y, 0, lc_f, 8)
        gmsh.model.geo.addPoint(cy_center_x, cy_center_y - cy_radious, 0, lc_f, 9)

        gmsh.model.geo.addCircleArc(6,5,7,5)
        gmsh.model.geo.addCircleArc(7,5,8,6)
        gmsh.model.geo.addCircleArc(8,5,9,7)
        gmsh.model.geo.addCircleArc(9,5,6,8)

        gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 12) #the circle
        gmsh.model.geo.synchronize()

        gmsh.model.geo.addPhysicalGroup(1, [5, 6, 7, 8], 12 )
        gmsh.model.setPhysicalName(1, 12, "circle")
        gmsh.model.geo.synchronize()

        # make the surface

        gmsh.model.geo.addPlaneSurface([10,12], 100) #the surface
        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(2, [100], 101)
        gmsh.model.setPhysicalName(2, 101, "surface")
        gmsh.model.geo.synchronize()

        gmsh.model.mesh.generate(2)
        gmsh.write("models/$name.msh")
        gmsh.finalize()
        model = GmshDiscreteModel("models/$name.msh")
    end

    if grid_type == "circle_circle" # circle - circle
    
        
        cy_center_x, cy_center_y, cy_inner_radious, cy_outer_radious, lc = p
        
        lc_f = lc
        lc = lc/4
        gmsh.model.add("$name")
        #lc = 1e-2
        #lc = lc/4
        gmsh.model.geo.addPoint(cy_center_x, cy_center_y, 0, lc, 1)
        gmsh.model.geo.addPoint(cy_center_x + cy_inner_radious, cy_center_y, 0, lc, 2)
        gmsh.model.geo.addPoint(cy_center_x , cy_center_y + cy_inner_radious, 0, lc, 3)
        gmsh.model.geo.addPoint(cy_center_x - cy_inner_radious, cy_center_y, 0, lc, 4)
        gmsh.model.geo.addPoint(cy_center_x, cy_center_y - cy_inner_radious, 0, lc, 5)


        # make the outer circle 
        gmsh.model.geo.addCircleArc(2,1,3,1)
        gmsh.model.geo.addCircleArc(3,1,4,2)
        gmsh.model.geo.addCircleArc(4,1,5,3)
        gmsh.model.geo.addCircleArc(5,1,2,4)


        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 10) #the outer circle
        gmsh.model.geo.synchronize()

        #gmsh.model.geo.addPhysicalGroup(1, [10], 11 )
        gmsh.model.geo.addPhysicalGroup(1, [1, 2, 3, 4], 11 )
        gmsh.model.setPhysicalName(1, 11, "inner")
        gmsh.model.geo.synchronize()

        # add the inner circle

        #lc_f = lc*4;
        gmsh.model.geo.addPoint(cy_center_x, cy_center_y, 0, lc_f, 6)
        gmsh.model.geo.addPoint(cy_center_x + cy_outer_radious, cy_center_y, 0, lc_f, 7)
        gmsh.model.geo.addPoint(cy_center_x , cy_center_y + cy_outer_radious, 0, lc_f, 8)
        gmsh.model.geo.addPoint(cy_center_x - cy_outer_radious, cy_center_y, 0, lc_f, 9)
        gmsh.model.geo.addPoint(cy_center_x, cy_center_y - cy_outer_radious, 0, lc_f, 10)

        gmsh.model.geo.addCircleArc(7,6,8,5)
        gmsh.model.geo.addCircleArc(8,6,9,6)
        gmsh.model.geo.addCircleArc(9,6,10,7)
        gmsh.model.geo.addCircleArc(10,6,7,8)

        gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 12) #the circle
        gmsh.model.geo.synchronize()

        gmsh.model.geo.addPhysicalGroup(1, [5, 6, 7, 8], 12 )
        gmsh.model.setPhysicalName(1, 12, "outer")
        gmsh.model.geo.synchronize()

        # make the surface

        gmsh.model.geo.addPlaneSurface([10,12], 100) #the surface
        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(2, [100], 101)
        gmsh.model.setPhysicalName(2, 101, "surface")
        gmsh.model.geo.synchronize()

        gmsh.model.mesh.generate(2)
        gmsh.write("models/$name.msh")
        gmsh.finalize()
        model = GmshDiscreteModel("models/$name.msh")
    end
    
    
    #======================= Nuevo dominio =======================#
    
    if grid_type == "rectangle_hole_square" # square - circle
        
        name, side_x, side_y, circ_center_x, circ_center_y, circ_radius, rec_base, rec_top, rec_left, rec_right, lc = p
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(name)
        #lc = 1e-2
        #lc = 1e-1
        gmsh.model.geo.addPoint(0, 0, 0, lc, 1)             # 1 vertice inferior izquierdo
        gmsh.model.geo.addPoint(side_x, 0,  0, lc, 2)       # 2 vértice inferior derecho
        gmsh.model.geo.addPoint(side_x, side_y, 0, lc, 3)   # 3 vértice superior derecho
        gmsh.model.geo.addPoint(0, side_y, 0, lc, 4)        # 4 vértice superior izquierdo

        # make the square boundary
        gmsh.model.geo.addLine(1, 2, 1) # 1 linea inferior
        gmsh.model.geo.addLine(2, 3, 2) # 2 línea lateral derecha
        gmsh.model.geo.addLine(3, 4, 3) # 3 linea superior
        gmsh.model.geo.addLine(4, 1, 4) # 4 linea lateral izquierda

        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 100) # the rectangle
        gmsh.model.geo.synchronize()

        #gmsh.model.geo.addPhysicalGroup(1, [10], 11 )
        gmsh.model.geo.addPhysicalGroup(1, [1, 2, 3, 4], 101 ) # grupo formado por cuatro lineas
        gmsh.model.setPhysicalName(1, 101, "ext")
        gmsh.model.geo.synchronize()

        # add the circle
        
        lc_f = lc/4;
        gmsh.model.geo.addPoint(circ_center_x, circ_center_y, 0, lc_f, 5)                   # 5 punto central en x,y
        gmsh.model.geo.addPoint(circ_center_x + circ_radius, circ_center_y, 0, lc_f, 6)     # 6 punto desplazado R en x & 0 en y
        gmsh.model.geo.addPoint(circ_center_x , circ_center_y + circ_radius, 0, lc_f, 7)    # 7 punto desplazado 0 en x & R en y
        gmsh.model.geo.addPoint(circ_center_x - circ_radius, circ_center_y, 0, lc_f, 8)     # 8 punto desplazado -R en x & 0 en y
        gmsh.model.geo.addPoint(circ_center_x, circ_center_y - circ_radius, 0, lc_f, 9)     # 9 punto desplazado 0 en x & -R en y

        gmsh.model.geo.addCircleArc(6,5,7,5) # 5 arco que une puntos 6,5,8
        gmsh.model.geo.addCircleArc(7,5,8,6) # 6 arco que une puntos 7,5,8
        gmsh.model.geo.addCircleArc(8,5,9,7) # 7 arco que une puntos 8,5,9
        gmsh.model.geo.addCircleArc(9,5,6,8) # 8 arco qu eune puntos 9,5,6

        gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 102) # the circle
        gmsh.model.geo.synchronize()

        gmsh.model.geo.addPhysicalGroup(1, [5, 6, 7, 8], 103 ) # grupo 1D formado por arcos
        gmsh.model.setPhysicalName(1, 103, "inner_circle")
        gmsh.model.geo.synchronize()
        
        # add the square
        
        lc_f = lc/4;
        gmsh.model.geo.addPoint(rec_left, rec_base, 0, lc_f, 10)    # 10 punto inferior izquierdo del rectangulo
        gmsh.model.geo.addPoint(rec_left, rec_top,  0, lc_f, 11)    # 11 punto superior izquierdo del rectangulo
        gmsh.model.geo.addPoint(rec_right, rec_top, 0, lc_f, 12)    # 12 punto superior derecho del rectangulo
        gmsh.model.geo.addPoint(rec_right, rec_base, 0, lc_f, 13)   # 13 punto inferior derecho del rectangulo

        # make the square boundary
        gmsh.model.geo.addLine(10, 11, 10)  # 10 linea lateral izquierda del rectangulo
        gmsh.model.geo.addLine(11, 12, 11)  # 11 linea superior del rectangulo
        gmsh.model.geo.addLine(12, 13, 12)  # 12 linea lateral derecha del rectangulo
        gmsh.model.geo.addLine(13, 10, 13)  # 13 linea inferior del rectangulo

        gmsh.model.geo.addCurveLoop([10, 11, 12, 13], 104) # the rectangle
        gmsh.model.geo.synchronize()

        #gmsh.model.geo.addPhysicalGroup(1, [10], 11 )
        gmsh.model.geo.addPhysicalGroup(1, [10, 11, 12, 13], 105 )
        gmsh.model.setPhysicalName(1, 105, "inner_square")
        gmsh.model.geo.synchronize()

        # make the surface

        gmsh.model.geo.addPlaneSurface([100, 102, 104], 1000) #the surface
        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(2, [1000], 1001)
        gmsh.model.setPhysicalName(2, 1001, "surface")
        gmsh.model.geo.synchronize()

        gmsh.model.mesh.generate(2)
        gmsh.write("models/$name.msh")
        gmsh.finalize()
        model = GmshDiscreteModel("models/$name.msh")
    end


    if grid_type == "rectangle_hole_square_v2" # square - circle
        
        name, side_x, side_y, circ_center_x, circ_center_y, circ_radius, rec_base, rec_top, rec_left, rec_right, lc, quad_state= p
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(name)
        #lc = 1e-2
        #lc = 1e-1
        gmsh.model.geo.addPoint(0, 0, 0, lc, 1)             # 1 vertice inferior izquierdo
        gmsh.model.geo.addPoint(side_x, 0,  0, lc, 2)       # 2 vértice inferior derecho
        gmsh.model.geo.addPoint(side_x, side_y, 0, lc, 3)   # 3 vértice superior derecho
        gmsh.model.geo.addPoint(0, side_y, 0, lc, 4)        # 4 vértice superior izquierdo

        # make the square boundary
        gmsh.model.geo.addLine(1, 2, 1) # 1 linea inferior
        gmsh.model.geo.addLine(2, 3, 2) # 2 línea lateral derecha
        gmsh.model.geo.addLine(3, 4, 3) # 3 linea superior
        gmsh.model.geo.addLine(4, 1, 4) # 4 linea lateral izquierda

        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 100) # the rectangle
        gmsh.model.geo.synchronize()

        #gmsh.model.geo.addPhysicalGroup(1, [10], 11 )
        gmsh.model.geo.addPhysicalGroup(1, [1, 2, 3, 4], 101 ) # grupo formado por cuatro lineas
        gmsh.model.setPhysicalName(1, 101, "ext")
        gmsh.model.geo.synchronize()

        # add the circle
        
        lc_f = lc/4;
        gmsh.model.geo.addPoint(circ_center_x, circ_center_y, 0, lc_f, 5)                   # 5 punto central en x,y
        gmsh.model.geo.addPoint(circ_center_x + circ_radius, circ_center_y, 0, lc_f, 6)     # 6 punto desplazado R en x & 0 en y
        gmsh.model.geo.addPoint(circ_center_x , circ_center_y + circ_radius, 0, lc_f, 7)    # 7 punto desplazado 0 en x & R en y
        gmsh.model.geo.addPoint(circ_center_x - circ_radius, circ_center_y, 0, lc_f, 8)     # 8 punto desplazado -R en x & 0 en y
        gmsh.model.geo.addPoint(circ_center_x, circ_center_y - circ_radius, 0, lc_f, 9)     # 9 punto desplazado 0 en x & -R en y

        gmsh.model.geo.addCircleArc(6,5,7,5) # 5 arco que une puntos 6,5,8
        gmsh.model.geo.addCircleArc(7,5,8,6) # 6 arco que une puntos 7,5,8
        gmsh.model.geo.addCircleArc(8,5,9,7) # 7 arco que une puntos 8,5,9
        gmsh.model.geo.addCircleArc(9,5,6,8) # 8 arco qu eune puntos 9,5,6

        gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 102) # the circle
        gmsh.model.geo.synchronize()

        gmsh.model.geo.addPhysicalGroup(1, [5, 6, 7, 8], 103 ) # grupo 1D formado por arcos
        gmsh.model.setPhysicalName(1, 103, "inner_circle")
        gmsh.model.geo.synchronize()
        
        # add the square
        
        lc_f = lc/4;
        gmsh.model.geo.addPoint(rec_left, rec_base, 0, lc_f, 10)    # 10 punto inferior izquierdo del rectangulo
        gmsh.model.geo.addPoint(rec_left, rec_top,  0, lc_f, 11)    # 11 punto superior izquierdo del rectangulo
        gmsh.model.geo.addPoint(rec_right, rec_top, 0, lc_f, 12)    # 12 punto superior derecho del rectangulo
        gmsh.model.geo.addPoint(rec_right, rec_base, 0, lc_f, 13)   # 13 punto inferior derecho del rectangulo

        # make the square boundary
        gmsh.model.geo.addLine(10, 11, 10)  # 10 linea lateral izquierda del rectangulo
        gmsh.model.geo.addLine(11, 12, 11)  # 11 linea superior del rectangulo
        gmsh.model.geo.addLine(12, 13, 12)  # 12 linea lateral derecha del rectangulo
        gmsh.model.geo.addLine(13, 10, 13)  # 13 linea inferior del rectangulo

        gmsh.model.geo.addCurveLoop([10, 11, 12, 13], 104) # the rectangle
        gmsh.model.geo.synchronize()

        #= +++++++++++++++++++++++++++++++++++++++++++++
            agregamos el siguiente bloque de código para incorporar vértices
            del rectangulo externo, circulo interno y rectangulo interno
           +++++++++++++++++++++++++++++++++++++++++++++ =#
        # gmsh.model.geo.addPhysicalGroup(0, [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13], 201 )
        gmsh.model.geo.addPhysicalGroup(0, [1, 2, 3, 4], 201 ) # 201 grupo formado por 4 puntos vértices del rectangulo externo
        gmsh.model.setPhysicalName(0, 201, "ext_Vertices")     # asignamos nombre al grupo
        gmsh.model.geo.synchronize()

        # gmsh.model.geo.addPhysicalGroup(0, [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13], 201 )
        gmsh.model.geo.addPhysicalGroup(0, [6, 7, 8, 9], 202 ) # 202 grupo formado por 4 puntos cartesianos del circulo interno 
        gmsh.model.setPhysicalName(0, 202, "intCirc_Vertices")
        gmsh.model.geo.synchronize()

        # gmsh.model.geo.addPhysicalGroup(0, [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13], 201 )
        gmsh.model.geo.addPhysicalGroup(0, [10, 11, 12, 13], 203 ) # 203 grupo formado por 4 vertices del rectangulo interno
        gmsh.model.setPhysicalName(0, 203, "intSquare_Vertices")
        gmsh.model.geo.synchronize()
        #= ++++++++++++++++++++++++++++++++++++++++++++ =#

        #gmsh.model.geo.addPhysicalGroup(1, [10], 11 )
        gmsh.model.geo.addPhysicalGroup(1, [10, 11, 12, 13], 105 )
        gmsh.model.setPhysicalName(1, 105, "inner_square")
        gmsh.model.geo.synchronize()

        # make the surface

        gmsh.model.geo.addPlaneSurface([100, 102, 104], 1000) #the surface
        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(2, [1000], 1001)
        gmsh.model.setPhysicalName(2, 1001, "surface")
        gmsh.model.geo.synchronize()

        if (quad_state=="true")
            gmsh.model.mesh.setRecombine(2, 1000) # for 2D quadrilaterals
        end

        gmsh.model.mesh.generate(2)
        gmsh.write("models/$name.msh")
        gmsh.finalize()
        model = GmshDiscreteModel("models/$name.msh")
    end

    #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Bloque de código agregado para crear un nuevo test que permita
        utilizar elementos cuadriláteros y mallas estructuradas.
       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

    if grid_type == "rectangle_hole_square_v3" # square - circle including vertices in Physical groups
        gmsh.option.setNumber("General.Terminal", 0)
        name, side_x, side_y, circ_center_x, circ_center_y, circ_radius, rec_base, rec_top, rec_left, rec_right, lc, quad_state = p
        gmsh.model.add(name)
        #lc = 1e-2
        #lc = 1e-1
        gmsh.model.geo.addPoint(0, 0, 0, lc, 1)             # 1 vertice inferior izquierdo
        gmsh.model.geo.addPoint(side_x, 0,  0, lc, 2)       # 2 vértice inferior derecho
        gmsh.model.geo.addPoint(side_x, side_y, 0, lc, 3)   # 3 vértice superior derecho
        gmsh.model.geo.addPoint(0, side_y, 0, lc, 4)        # 4 vértice superior izquierdo

        # make the square boundary

        gmsh.model.geo.addLine(1, 2, 1) # 1 linea inferior
        gmsh.model.geo.addLine(2, 3, 2) # 2 línea lateral derecha
        gmsh.model.geo.addLine(3, 4, 3) # 3 linea superior
        gmsh.model.geo.addLine(4, 1, 4) # 4 linea lateral izquierda

        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 100) # the rectangle
        gmsh.model.geo.synchronize()

        #gmsh.model.geo.addPhysicalGroup(1, [10], 11 )
        
        gmsh.model.geo.addPhysicalGroup(1, [1, 2, 3, 4], 101 ) # grupo formado por cuatro lineas
        gmsh.model.setPhysicalName(1, 101, "ext")
        gmsh.model.geo.synchronize()

        lc_f = lc/4;

        # add the circle points
        gmsh.model.geo.addPoint(circ_center_x, circ_center_y, 0, lc_f, 5)                   # 5 punto central en x,y
        gmsh.model.geo.addPoint(circ_center_x + circ_radius, circ_center_y, 0, lc_f, 6)     # 6 punto desplazado R en x & 0 en y
        gmsh.model.geo.addPoint(circ_center_x , circ_center_y + circ_radius, 0, lc_f, 7)    # 7 punto desplazado 0 en x & R en y
        gmsh.model.geo.addPoint(circ_center_x - circ_radius, circ_center_y, 0, lc_f, 8)     # 8 punto desplazado -R en x & 0 en y
        gmsh.model.geo.addPoint(circ_center_x, circ_center_y - circ_radius, 0, lc_f, 9)     # 9 punto desplazado 0 en x & -R en y

        # add the square points
        gmsh.model.geo.addPoint(rec_left, rec_base, 0, lc_f, 10)    # 10 punto inferior izquierdo del rectangulo
        gmsh.model.geo.addPoint(rec_left, rec_top,  0, lc_f, 11)    # 11 punto superior izquierdo del rectangulo
        gmsh.model.geo.addPoint(rec_right, rec_top, 0, lc_f, 12)    # 12 punto superior derecho del rectangulo
        gmsh.model.geo.addPoint(rec_right, rec_base, 0, lc_f, 13)   # 13 punto inferior derecho del rectangulo

        #= +++++++++++++++++++++++++++++++++++++++++++++
            agregamos el siguiente bloque de código para incorporar vértices
            del rectangulo externo, circulo interno y rectangulo interno
            a un nuevo grupo para aplicar condiciones de borde
           +++++++++++++++++++++++++++++++++++++++++++++ =#
        gmsh.model.geo.addPhysicalGroup(0, [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13], 201 )
        gmsh.model.setPhysicalName(0, 201, "ext_Vertices")
        gmsh.model.geo.synchronize()
        # ===============================================

        # add the circle
        gmsh.model.geo.addCircleArc(6,5,7,5) # 5 arco que une puntos 6,5,8
        gmsh.model.geo.addCircleArc(7,5,8,6) # 6 arco que une puntos 7,5,8
        gmsh.model.geo.addCircleArc(8,5,9,7) # 7 arco que une puntos 8,5,9
        gmsh.model.geo.addCircleArc(9,5,6,8) # 8 arco qu eune puntos 9,5,6

        gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 102) # the circle
        gmsh.model.geo.synchronize()

        gmsh.model.geo.addPhysicalGroup(1, [5, 6, 7, 8], 103 )
        gmsh.model.setPhysicalName(1, 103, "inner_circle")
        gmsh.model.geo.synchronize()

        # add the square

        # make the square boundary
        gmsh.model.geo.addLine(10, 11, 10)  # 10 linea lateral izquierda del rectangulo
        gmsh.model.geo.addLine(11, 12, 11)  # 11 linea superior del rectangulo
        gmsh.model.geo.addLine(12, 13, 12)  # 12 linea lateral derecha del rectangulo
        gmsh.model.geo.addLine(13, 10, 13)  # 13 linea inferior del rectangulo

        gmsh.model.geo.addCurveLoop([10, 11, 12, 13], 104) # the rectangle
        gmsh.model.geo.synchronize()

        gmsh.model.geo.addPhysicalGroup(1, [10, 11, 12, 13], 105 )
        gmsh.model.setPhysicalName(1, 105, "inner_square")
        gmsh.model.geo.synchronize()

        # make the surface

        gmsh.model.geo.addPlaneSurface([100, 102, 104], 1000) # the surface
        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(2, [1000], 1001)
        gmsh.model.setPhysicalName(2, 1001, "surface")
        gmsh.model.geo.synchronize()

        if (quad_state=="true")
            gmsh.model.mesh.setRecombine(2, 1000) # for 2D quadrilaterals
        end

        gmsh.model.mesh.generate(2)
        gmsh.write("models/$name.msh")
        gmsh.finalize()
        model = GmshDiscreteModel("models/$name.msh")
    end


    #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#

    
    if grid_type == "rectangle_point"  # simple square with a point at center to define Dirac
        
        name, side_x, side_y, lc = p
        #first we build the rectangular boundary: 
        #gmsh.initialize()
        #gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("$name")
        #lc = 1e-2
        #lc = 1e-1
        lc_x = lc
        lc_y = lc#*side_y/side_x
        gmsh.model.geo.addPoint(0, 0, 0, lc_x, 1)
        gmsh.model.geo.addPoint(side_x, 0,  0, lc_x, 2)
        gmsh.model.geo.addPoint(side_x, side_y, 0, lc_y, 3)
        gmsh.model.geo.addPoint(0, side_y, 0, lc_y, 4)
        gmsh.model.geo.addPoint(side_x/2, side_y/2, 0, lc/10, 5)
        gmsh.model.geo.synchronize()
        
        gmsh.model.geo.addPhysicalGroup(0, [5], 1 )
        gmsh.model.setPhysicalName(0, 1, "point")
        gmsh.model.geo.synchronize()

        # make the square boundary
        gmsh.model.geo.addLine(1, 2, 1)
        gmsh.model.geo.addLine(2, 3, 2)
        gmsh.model.geo.addLine(3, 4, 3)
        gmsh.model.geo.addLine(4, 1, 4)

        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 10) #the rectangle
        gmsh.model.geo.synchronize()

        #gmsh.model.geo.addPhysicalGroup(1, [10], 11 )
        gmsh.model.geo.addPhysicalGroup(1, [1, 2, 3, 4], 11 )
        gmsh.model.setPhysicalName(1, 11, "ext")
        gmsh.model.geo.synchronize()

        # make the surface

        gmsh.model.geo.addPlaneSurface([10], 100) #the surface
        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(2, [100], 101)
        gmsh.model.setPhysicalName(2, 101, "surface")
        gmsh.model.geo.synchronize()

        gmsh.model.mesh.generate(2)
        gmsh.write("models/$name.msh")
        gmsh.finalize()
        model = GmshDiscreteModel("models/$name.msh")
 
    end



    #gmsh.finalize()
    return model
end