#include "mfem.hpp"
#include "TLIntegrators.hpp"
#include "TLStressStrain.hpp"

// Inclusion radius and radial displacement 
double r_inc, u_r;

void u_inc(const mfem::Vector& x, double factor, mfem::Vector& u) 
{
    u.SetSize(2);
    u[0] = u_r * x[0] / r_inc;
    u[1] = u_r * x[1] / r_inc;
    u *= factor;
}

int main(int argc, char** argv)
{
    // Material parameters
    double a, A1, A2, A3, A4, A5, A6;
    
    // I/O parameters
    std::string MeshFile, ResultFile;

    mfem::OptionsParser args(argc, argv);
    args.AddOption(&MeshFile, "-mf", "--MeshFile", "Mesh File");
    args.AddOption(&ResultFile, "-rf", "--ResultFile", "Result File");
    args.AddOption(&a, "-a", "--a", "a");
    args.AddOption(&A1, "-A1", "--A1", "A1");
    args.AddOption(&A2, "-A2", "--A2", "A2");
    args.AddOption(&A3, "-A3", "--A3", "A3");
    args.AddOption(&A4, "-A4", "--A4", "A4");
    args.AddOption(&A5, "-A5", "--A5", "A5");
    args.AddOption(&A6, "-A6", "--A6", "A6");
    args.AddOption(&r_inc, "-r_inc", "--r_inc", "Inclusion radius");
    args.AddOption(&u_r, "-u_r", "--u_r", "Inclusion radial displacement");
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

    auto a_coeff = mfem::ConstantCoefficient(a);
    auto A1_coeff = mfem::ConstantCoefficient(A1);
    auto A2_coeff = mfem::ConstantCoefficient(A2);
    auto A3_coeff = mfem::ConstantCoefficient(A3);
    auto A4_coeff = mfem::ConstantCoefficient(A4);
    auto A5_coeff = mfem::ConstantCoefficient(A5);
    auto A6_coeff = mfem::ConstantCoefficient(A6);

    mfem::Array<int> ess_tdofs, tmp_tdofs;
    mfem::Array<int> inclusion({1, 0});
    mfem::Array<int> far_field({0, 1});

    u_space.GetEssentialTrueDofs(far_field, tmp_tdofs);
    ess_tdofs.Append(tmp_tdofs);
    
    u_space.GetEssentialTrueDofs(inclusion, tmp_tdofs);
    ess_tdofs.Append(tmp_tdofs);

    auto inclusion_bc_coeff = mfem::VectorFunctionCoefficient(dim, u_inc);

    auto B = mfem::NonlinearForm(&u_space);
    B.AddDomainIntegrator(new FungExponentialIntegrator(a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                                                        A4_coeff, A5_coeff, A6_coeff));
    B.SetEssentialTrueDofs(ess_tdofs);

    auto prec = mfem::UMFPackSolver();
    auto ns = mfem::NewtonSolver();
    ns.SetOperator(B);
    ns.SetPreconditioner(prec);
    ns.SetRelTol(1e-14);
    ns.SetAbsTol(1e-8);
    ns.SetMaxIter(40);
    ns.SetPrintLevel(1);

    int N_increments = 100;
    for (int i=0; i<N_increments; i++)
    {
        mfem::out << "Solving increment " << (i+1) << " out of " << N_increments << " \n";
        // Use "time" to incrementally apply bc
        inclusion_bc_coeff.SetTime(static_cast<double>(i+1)/N_increments); 
        u.ProjectBdrCoefficient(inclusion_bc_coeff, inclusion);
        ns.Mult(f, u);
    }

    std::ofstream file(ResultFile);
    file.precision(16);
    mesh.PrintVTK(file, 0);
    u.SaveVTK(file, "u", 0);
    file.close();

    return 0;
}
