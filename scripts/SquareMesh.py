import gmsh
from math import sqrt

SideLength = 1.
NumberOfElements = 4096
tol = 1e-6 # Do not make tolerance smaller

def main():

    gmsh.initialize()

    # Create square and two lines that intersect at its center
    square = gmsh.model.occ.addRectangle(0., 0., 0., SideLength, SideLength)
    p_left   = gmsh.model.occ.addPoint(0.0, 0.5, 0.0)
    p_right  = gmsh.model.occ.addPoint(1.0, 0.5, 0.0)
    p_bottom = gmsh.model.occ.addPoint(0.5, 0.0, 0.0)
    p_top    = gmsh.model.occ.addPoint(0.5, 1.0, 0.0)
    hline = gmsh.model.occ.addLine(p_left, p_right)
    vline = gmsh.model.occ.addLine(p_bottom, p_top)
    gmsh.model.occ.fragment([(2, square)], [(1, hline), (1, vline)])
    gmsh.model.occ.synchronize()

    surfaces = gmsh.model.getEntities(2)
    edges = gmsh.model.getEntities(1)

    left, right, bottom, top, hsplit, vsplit = [], [], [], [], [], []

    for dim, tag in edges:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)

        if abs(xmin) < tol and abs(xmax) < tol:
            left.append(tag)
        elif abs(xmin - 1.0) < tol and abs(xmax - 1.0) < tol:
            right.append(tag)
        elif abs(ymin) < tol and abs(ymax) < tol:
            bottom.append(tag)
        elif abs(ymin - 1.0) < tol and abs(ymax - 1.0) < tol:
            top.append(tag)
        elif abs(ymin - 0.5) < tol and abs(ymax - 0.5) < tol:
            hsplit.append(tag)
        elif abs(xmin - 0.5) < tol and abs(xmax - 0.5) < tol:
            vsplit.append(tag)

    gmsh.model.addPhysicalGroup(1, bottom, tag=1, name="Bottom")
    gmsh.model.addPhysicalGroup(1, right, tag=2, name="Right")
    gmsh.model.addPhysicalGroup(1, top, tag=3, name="Top")
    gmsh.model.addPhysicalGroup(1, left, tag=4, name="Left")
    gmsh.model.addPhysicalGroup(1, hsplit, tag=5, name="hSplit")
    gmsh.model.addPhysicalGroup(1, vsplit, tag=6, name="vSplit")

    domain = [s[1] for s in surfaces]
    gmsh.model.addPhysicalGroup(2, domain, tag=8, name="Domain")

    # Set mesh size at all mesh points
    ElementsOnEdge = sqrt(NumberOfElements)
    MeshSize = SideLength/ElementsOnEdge
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), MeshSize)

    # Make the mesh quadrilateral
    for dim, tag in surfaces:
        gmsh.model.mesh.setTransfiniteSurface(tag)
    
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.model.mesh.generate(dim=2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write("../meshes/Square.msh")
    gmsh.clear()

    gmsh.finalize()

if __name__ == "__main__":
    main()
