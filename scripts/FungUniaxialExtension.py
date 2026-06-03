import subprocess

A0 = 10.
a = 200e3
strains = [0.01*i for i in range(1,31)]

def main():

    CommandLineInput = ["../build/problems/FUE"]
    CommandLineInput.extend(["--MeshFile", "../meshes/Square.msh"])
    CommandLineInput.extend(["--ResultFile", " "])
    CommandLineInput.extend(["-a", str(a)])
    CommandLineInput.extend(["-A1", str(A0)])
    CommandLineInput.extend(["-A2", str(A0)])
    CommandLineInput.extend(["-A3", str(A0/2.)])
    CommandLineInput.extend(["-A4", str(A0/2.)])
    CommandLineInput.extend(["-A5", "0."])
    CommandLineInput.extend(["-A6", "0."])
    CommandLineInput.extend(["-u_x", str(strains[0])])

    processes = []
    for i in range(0, len(strains)):
        CommandLineInput[4] = str("../results/Fung_Uniaxial_Extension_Result_" + str(i) + ".vtk")
        CommandLineInput[20] = str(strains[i])
        process = subprocess.Popen(CommandLineInput)
        processes.append(process)

    for process in processes:
        returncode = process.wait()
        assert int(returncode) == 0, "\033[31mProcess Failed\033[0m"

if __name__ == "__main__":
    main()
