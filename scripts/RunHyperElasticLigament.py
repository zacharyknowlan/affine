import subprocess
from LigamentMesh import l

Lambda = 0.4
mu = 0.2

strains = [0.05*i for i in range(1,9)]
displacements = [l*strain for strain in strains]

def main():

    CommandLineInput = ["../build/problems/HEL"]
    CommandLineInput.extend(["--MeshFile", "../meshes/HalfLigament.msh"])
    CommandLineInput.extend(["--ResultFile", " "])
    CommandLineInput.extend(["--lambda", str(Lambda)])
    CommandLineInput.extend(["--mu", str(mu)])
    CommandLineInput.extend(["--u", " "])

    print("Running Hyper Elastic Ligament Cases...")
    for i in range(0, len(displacements)):
        CommandLineInput[4] = str("../results/HyperElasticLigament_" + str(i) + ".vtk")
        CommandLineInput[10] = str(displacements[i])
        _ = subprocess.run(CommandLineInput)

if __name__ == "__main__":
    main()
