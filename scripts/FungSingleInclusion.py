import subprocess
from SingleInclusionMesh import InclusionRadius

def main():

    E = 10.
    
    CommandLineInput = ["../build/example_problems/FungSingleInclusion"]
    CommandLineInput.extend(["--MeshFile", "../meshes/SingleInclusionFarField.msh"])
    CommandLineInput.extend(["--ResultFile", "../results/FungSingleInclusion.vtk"])
    CommandLineInput.extend(["-a", str(5.*E)])
    CommandLineInput.extend(["-A1", str(E)])
    CommandLineInput.extend(["-A2", str(E)])
    CommandLineInput.extend(["-A3", str(E/2.)])
    CommandLineInput.extend(["-A4", str(E/2.)])
    CommandLineInput.extend(["-A5", "0."])
    CommandLineInput.extend(["-A6", "0."])
    CommandLineInput.extend(["-r_inc", str(InclusionRadius)])
    CommandLineInput.extend(["-u_r", str(-0.32*InclusionRadius)])

    _ = subprocess.run(CommandLineInput)

if __name__ == "__main__":
    main()
