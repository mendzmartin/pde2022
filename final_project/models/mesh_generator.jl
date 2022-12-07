function make_model(grid_type,p)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    # one dimensional grid

    if grid_type == "1D"
        x0, x1, lc = p
        
        gmsh.model.add("1D")
        
        gmsh.model.geo.addPoint(x0, 0, 0, lc, 1)
        gmsh.model.geo.addPoint(x1, 0, 0, lc, 2)

        gmsh.model.geo.addLine(1, 2, 1)
        gmsh.model.geo.synchronize()
        
        gmsh.model.geo.addPhysicalGroup(0, [1], 3 )
        gmsh.model.setPhysicalName(0, 3, "left")
        gmsh.model.geo.addPhysicalGroup(0, [2], 4 )
        gmsh.model.setPhysicalName(0, 4, "right")
        
        gmsh.model.geo.addPhysicalGroup(1, [1], 1 )
        gmsh.model.setPhysicalName(1, 1, "segment")
        gmsh.model.geo.synchronize()


        gmsh.model.mesh.generate(1)
        gmsh.write("../models/1D.msh")
        
        gmsh.finalize()
        model = GmshDiscreteModel("../models/1D.msh")
    end

end