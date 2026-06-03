import subprocess
from RunHyperElasticLigament import Lambda, mu, displacements

def main():

    c = 20.

    CommandLineInput = ["../build/problems/FL"]
    CommandLineInput.extend(["--MeshFile", "../meshes/HalfLigament.msh"])
    CommandLineInput.extend(["--ResultFile", " "])
    CommandLineInput.extend(["--a", str(1./c)])
    CommandLineInput.extend(["--A1", str(c*(2.*mu+Lambda))])
    CommandLineInput.extend(["--A2", str(c*(2.*mu+Lambda))])
    CommandLineInput.extend(["--A3", str(c*Lambda)])
    CommandLineInput.extend(["--A4", str(c*2.*mu)])
    CommandLineInput.extend(["--A5", "0."])
    CommandLineInput.extend(["--A6", "0."])
    CommandLineInput.extend(["--u", " "])

    print("Running Fung Ligament Cases...")
    for i in range(0, len(displacements)):
        CommandLineInput[4] = str("../results/FungLigament_" + str(i) + ".vtk")
        CommandLineInput[20] = str(displacements[i])
        _ = subprocess.run(CommandLineInput)

if __name__ == "__main__":
    main()
