#include "mfem.hpp"
#include "affine.hpp"
#include <cmath>

double p_mag; // Internal pressure of the artery

void p_int(const mfem::Vector& x, double pseudo_time, mfem::Vector& T) 
{
    T.SetSize(2);
    double r = std::sqrt(std::pow(x[0], 2) + std::pow(x[1], 2));
    T[0] = (x[0]/r)*p_mag;
    T[1] = (x[1]/r)*p_mag;
    T *= pseudo_time;
}

int main(int argc, char** argv)
{
    // Material parameters
    double lambda, mu;
    
    // I/O parameters
    std::string MeshFile, ResultFile;

    mfem::OptionsParser args(argc, argv);
    args.AddOption(&MeshFile, "-mf", "--MeshFile", "Mesh File");
    args.AddOption(&ResultFile, "-rf", "--ResultFile", "Result File");
    args.AddOption(&lambda, "-lambda", "--lambda", "lambda");
    args.AddOption(&mu, "-mu", "--mu", "mu");
    args.AddOption(&p_mag, "-p", "--p", "Internal pressure");
    args.Parse();
    if (!args.Good())
    {
       args.PrintUsage(std::cout);
       return 1;
    }

    mfem::Mesh mesh = mfem::Mesh(MeshFile.c_str(), 1, 1);
    int dim = mesh.Dimension();

    auto u_ec = mfem::H1_FECollection(2, dim, mfem::BasisType::GaussLobatto); 
    auto u_space = mfem::FiniteElementSpace(&mesh, &u_ec, dim);

    mfem::GridFunction u(&u_space);
    u = 0.;
    mfem::GridFunction f(&u_space);
    f = 0.;

    auto lambda_coeff = mfem::ConstantCoefficient(lambda);
    auto mu_coeff = mfem::ConstantCoefficient(mu);

    mfem::Array<int> ess_tdofs, tmp_tdofs;
    mfem::Array<int> horizontal_edge({1, 0, 0, 0});
    mfem::Array<int> vertical_edge({0, 0, 1, 0});
    mfem::Array<int> inner_arc({0, 0, 0, 1});

    u_space.GetEssentialTrueDofs(horizontal_edge, tmp_tdofs, 1);
    ess_tdofs.Append(tmp_tdofs);
    
    u_space.GetEssentialTrueDofs(vertical_edge, tmp_tdofs, 0);
    ess_tdofs.Append(tmp_tdofs);

    auto T = mfem::VectorFunctionCoefficient(2, p_int);
    
    auto B = mfem::NonlinearForm(&u_space);
    B.AddDomainIntegrator(new HyperElasticIntegrator(mu_coeff, lambda_coeff));
    B.AddBoundaryIntegrator(new PK1TractionIntegrator(T), inner_arc);
    B.SetEssentialTrueDofs(ess_tdofs);

    auto prec = mfem::UMFPackSolver();
    auto ns = mfem::NewtonSolver();
    ns.SetOperator(B);
    ns.SetPreconditioner(prec);
    ns.SetRelTol(1e-12);
    ns.SetAbsTol(1e-8);
    ns.SetMaxIter(40);
    ns.SetPrintLevel(0);

    int N_increments = 50;
    for (int i=0; i<N_increments; i++)
    {
        // Pseudo time is fraction of applied load
        T.SetTime(static_cast<double>(i+1)/N_increments); 
        ns.Mult(f, u);
    }
    
    auto dg_ec = mfem::DG_FECollection(0, dim, mfem::BasisType::GaussLegendre);
    auto dg_tensor_space = mfem::FiniteElementSpace(&mesh, &dg_ec, dim*dim);
    auto dg_scalar_space = mfem::FiniteElementSpace(&mesh, &dg_ec, 1);

    auto E = mfem::GridFunction(&dg_tensor_space);
    CalcGreenLagrangeStrain(u, E);

    auto sigma = mfem::GridFunction(&dg_tensor_space);
    CalcHyperElasticCauchyStress(u, E, mu_coeff, lambda_coeff, sigma);

    auto sigma_VM = mfem::GridFunction(&dg_scalar_space);
    CalcVonMisesStress(sigma, sigma_VM);

    std::ofstream file(ResultFile);
    file.precision(16);
    mesh.PrintVTK(file, 0);
    u.SaveVTK(file, "u", 0);
    E.SaveVTK(file, "E", 0);
    sigma.SaveVTK(file, "sigma", 0);
    sigma_VM.SaveVTK(file, "sigma_VM", 0);
    file.close();

    return 0;
}
