import subprocess

Lambda = 0.4
mu = 0.2

pressures = [(0.016 + 0.001*i) for i in range(0,9)] # in MPa

def main():

    CommandLineInput = ["../build/problems/HEA"]
    CommandLineInput.extend(["--MeshFile", "../meshes/QuarterArtery.msh"])
    CommandLineInput.extend(["--ResultFile", " "])
    CommandLineInput.extend(["--lambda", str(Lambda)])
    CommandLineInput.extend(["--mu", str(mu)])
    CommandLineInput.extend(["--p", " "])

    print("Running Hyper Elastic Artery Cases...")
    for i in range(0, len(pressures)):
        CommandLineInput[4] = str("../results/HyperElasticArtery_" + str(i) + ".vtk")
        CommandLineInput[10] = str(pressures[i])
        _ = subprocess.run(CommandLineInput)

if __name__ == "__main__":
    main()
