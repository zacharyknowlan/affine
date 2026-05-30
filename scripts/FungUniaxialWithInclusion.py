import subprocess

E = 10.
strain = 0.01

def main():

    CommandLineInput = ["../build/example_problems/FungUniaxialExtension"]
    CommandLineInput.extend(["--MeshFile", "../meshes/UnitSquareWithInclusion.msh"])
    CommandLineInput.extend(["--ResultFile", "../results/SquareWithInclusion.vtk"])
    CommandLineInput.extend(["-a", str(5.*E)])
    CommandLineInput.extend(["-A1", str(E)])
    CommandLineInput.extend(["-A2", str(E)])
    CommandLineInput.extend(["-A3", str(E/2.)])
    CommandLineInput.extend(["-A4", str(E/2.)])
    CommandLineInput.extend(["-A5", "0."])
    CommandLineInput.extend(["-A6", "0."])
    CommandLineInput.extend(["-u_x", str(strain)])

    _ = subprocess.run(CommandLineInput)

if __name__ == "__main__":
    main()
