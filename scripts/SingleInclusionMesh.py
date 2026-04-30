import gmsh

InclusionRadius = 0.001
FarFieldRadius = 1e3 * InclusionRadius
MeshFile = "../meshes/SingleInclusionFarField.msh"
MeshSizes = [InclusionRadius/20., FarFieldRadius/30.]

def main():

    gmsh.initialize()

    Inclusion = gmsh.model.occ.addDisk(0., 0., 0., InclusionRadius, InclusionRadius)
    FarField = gmsh.model.occ.addDisk(0., 0., 0., FarFieldRadius, FarFieldRadius)
    Cut = gmsh.model.occ.cut([(2, FarField)], [(2, Inclusion)])
    gmsh.model.occ.synchronize()

    Boundary = gmsh.model.getBoundary(Cut[0], oriented=False)
    gmsh.model.mesh.setSize([(0, Boundary[0][1])], MeshSizes[0])
    gmsh.model.mesh.setSize([(0, Boundary[1][1])], MeshSizes[1])

    gmsh.model.addPhysicalGroup(1, [Boundary[0][1]], tag=1, name="Inclusion")
    gmsh.model.addPhysicalGroup(1, [Boundary[1][1]], tag=2, name="Far Field")
    Surface = [Entity[1] for Entity in gmsh.model.getEntities(dim=2)]
    gmsh.model.addPhysicalGroup(2, Surface, tag=3, name="Surface")

    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.model.mesh.generate(dim=2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(MeshFile)

    gmsh.finalize()

if __name__ == "__main__":
    main()
