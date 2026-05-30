import gmsh
import numpy as np

def main():

    gmsh.initialize()
    gmsh.model.add("UnitSquareWithInclusion")

    L = 1.
    cx, cy = L/2., L/2.
    r = 0.1

    square = gmsh.model.occ.addRectangle(0., 0., 0., L, L)
    circle = gmsh.model.occ.addDisk(cx, cy, 0., r, r)
    cut, _ = gmsh.model.occ.cut( [(2, square)], [(2, circle)], removeObject=True, removeTool=True)

    p_left = gmsh.model.occ.addPoint(0., 0.5, 0.)
    p_left_inc = gmsh.model.occ.addPoint(cx-r, 0.5, 0.)
    p_right_inc = gmsh.model.occ.addPoint(cx+r, 0.5, 0.)
    p_right  = gmsh.model.occ.addPoint(1., 0.5, 0.)
    line_left = gmsh.model.occ.addLine(p_left, p_left_inc)
    line_right = gmsh.model.occ.addLine(p_right_inc, p_right)
    gmsh.model.occ.fragment(cut, [(1, line_left), (1, line_right)])
    
    gmsh.model.occ.synchronize()

    edges = gmsh.model.getEntities(1)

    left, right, bottom, top, centerline, inc = [], [], [], [], [], []
    tol = 1e-6 # Do not make tolerance smaller

    for dim, tag in edges:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)

        com = gmsh.model.occ.getCenterOfMass(dim, tag)
        dist = np.sqrt((com[0] - cx)**2 + (com[1] - cy)**2)

        if abs(xmin) < tol and abs(xmax) < tol:
            left.append(tag)
        elif abs(xmin - 1.) < tol and abs(xmax - 1.) < tol:
            right.append(tag)
        elif abs(ymin) < tol and abs(ymax) < tol:
            bottom.append(tag)
        elif abs(ymin - 1.) < tol and abs(ymax - 1.) < tol:
            top.append(tag)
        elif abs(ymin - L/2.) < tol and abs(ymax - L/2.) < tol and (xmax > (cx + r) or xmin < (cx - r)):
            centerline.append(tag)
        elif abs(dist - r) < tol:
            inc.append(tag)

    gmsh.model.addPhysicalGroup(1, bottom, tag=1, name="Bottom")
    gmsh.model.addPhysicalGroup(1, right, tag=2, name="Right")
    gmsh.model.addPhysicalGroup(1, top, tag=3, name="Top")
    gmsh.model.addPhysicalGroup(1, left, tag=4, name="Left")
    gmsh.model.addPhysicalGroup(1, centerline, tag=5, name="Centerline")
    gmsh.model.addPhysicalGroup(1, inc, tag=6, name="Inclusion")

    surfaces = gmsh.model.getEntities(2)
    domain = [s[1] for s in surfaces]
    gmsh.model.addPhysicalGroup(2, domain, tag=7, name="Domain")

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.005)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.01)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(dim=2)
    gmsh.write("../meshes/UnitSquareWithInclusion.msh")
    
    gmsh.finalize()

if __name__ == "__main__":
    main()
