import subprocess
from SingleInclusionMesh import InclusionRadius

def main():

    A0 = 10.
    a = 200e3
    
    CommandLineInput = ["../build/problems/FSI"]
    CommandLineInput.extend(["--MeshFile", "../meshes/SingleInclusionFarField.msh"])
    CommandLineInput.extend(["--ResultFile", "../results/FungSingleInclusion.vtk"])
    CommandLineInput.extend(["-a", str(a)])
    CommandLineInput.extend(["-A1", str(A0)])
    CommandLineInput.extend(["-A2", str(A0)])
    CommandLineInput.extend(["-A3", str(A0/2.)])
    CommandLineInput.extend(["-A4", str(A0/2.)])
    CommandLineInput.extend(["-A5", "0."])
    CommandLineInput.extend(["-A6", "0."])
    CommandLineInput.extend(["-r_inc", str(InclusionRadius)])
    CommandLineInput.extend(["-u_r", str(-0.32*InclusionRadius)])

    _ = subprocess.run(CommandLineInput)

if __name__ == "__main__":
    main()
