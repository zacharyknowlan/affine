import subprocess
from RunHyperElasticArtery import Lambda, mu, pressures

def main():

    c = 20.

    CommandLineInput = ["../build/problems/FA"]
    CommandLineInput.extend(["--MeshFile", "../meshes/QuarterArtery.msh"])
    CommandLineInput.extend(["--ResultFile", " "])
    CommandLineInput.extend(["--a", str(1./c)])
    CommandLineInput.extend(["--A1", str(c*(2.*mu+Lambda))])
    CommandLineInput.extend(["--A2", str(c*(2.*mu+Lambda))])
    CommandLineInput.extend(["--A3", str(c*Lambda)])
    CommandLineInput.extend(["--A4", str(c*2.*mu)])
    CommandLineInput.extend(["--A5", "0."])
    CommandLineInput.extend(["--A6", "0."])
    CommandLineInput.extend(["--p", " "])

    print("Running Fung Artery Cases...")
    for i in range(0, len(pressures)):
        CommandLineInput[4] = str("../results/FungArtery_" + str(i) + ".vtk")
        CommandLineInput[20] = str(pressures[i])
        _ = subprocess.run(CommandLineInput)

if __name__ == "__main__":
    main()
