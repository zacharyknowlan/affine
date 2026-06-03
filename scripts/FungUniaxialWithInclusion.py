import subprocess

A0 = 10.
a = 200e3
strain = 0.01

def main():

    CommandLineInput = ["../build/problems/FUE"]
    CommandLineInput.extend(["--MeshFile", "../meshes/UnitSquareWithInclusion.msh"])
    CommandLineInput.extend(["--ResultFile", "../results/SquareWithInclusion.vtk"])
    CommandLineInput.extend(["-a", str(a)])
    CommandLineInput.extend(["-A1", str(A0)])
    CommandLineInput.extend(["-A2", str(A0)])
    CommandLineInput.extend(["-A3", str(A0/2.)])
    CommandLineInput.extend(["-A4", str(A0/2.)])
    CommandLineInput.extend(["-A5", "0."])
    CommandLineInput.extend(["-A6", "0."])
    CommandLineInput.extend(["-u_x", str(strain)])

    _ = subprocess.run(CommandLineInput)

if __name__ == "__main__":
    main()
