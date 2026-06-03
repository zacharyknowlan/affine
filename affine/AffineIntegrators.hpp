#ifndef AFFINE_INTEGRATORS_HPP
#define AFFINE_INTEGRATORS_HPP

#include "mfem.hpp"

class FungExponentialIntegrator : public mfem::NonlinearFormIntegrator
{
    private:

        // Matrices used for element residual/gradient computation
        mfem::DenseMatrix Gradu_, F_, S_, E_, A_, dSdu_;
        mfem::DenseMatrix dNdeta_, dNdx_, dLambdadE_, tmp1_, tmp2_;

        // Material properties and stiffening exponent
        double a_, A1_, A2_, A3_, A4_, A5_, A6_, Lambda_;

        // Material property coefficients
        mfem::Array<mfem::Coefficient*> mat_props_;

        // Get the Voigt index corresponding to i and j
        int MapVoigt2D(const int i, const int j);

        // Compute the Green-Lagrange strain tensor
        // Result is stored in E_
        void ComputeE();

        // Compute the energy exponent Lambda
        // Result is stored in Lambda_
        void ComputeLambda();

        // Compute the derivative of the energy exponential with 
        // respect to the Green-Lagrange strain
        // Result is stored in dLambdadE_
        void ComputedLambdadE();

        // Compute the constant portion of the Voigt stiffness tensor
        // The result is stored in A_
        void ComputeA();

    public:

        FungExponentialIntegrator(mfem::Coefficient& a, mfem::Coefficient& A1, mfem::Coefficient& A2,
                                    mfem::Coefficient& A3, mfem::Coefficient& A4, mfem::Coefficient& A5,
                                    mfem::Coefficient& A6);

        void AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                    const mfem::Vector& elfun, mfem::Vector& elvect) override;

        void AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                    const mfem::Vector& elfun, mfem::DenseMatrix& elmat) override;
};

class HyperElasticIntegrator : public mfem::NonlinearFormIntegrator
{
    private:

        // Matrices used for element residual/gradient computation
        mfem::DenseMatrix Gradu_, F_, S_, dSdu_, tmp1_, dNdeta_, dNdx_;

        // Material property values
        double mu_, lambda_;

        // Material property coefficients
        mfem::Array<mfem::Coefficient*> mat_props_;

    public:

        HyperElasticIntegrator(mfem::Coefficient& mu, mfem::Coefficient& lambda);

        void AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                    const mfem::Vector& elfun, mfem::Vector& elvect) override;

        void AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                    const mfem::Vector& elfun, mfem::DenseMatrix& elmat) override;

};

class PK1TractionIntegrator : public mfem::NonlinearFormIntegrator
{
    private:

        // Matrices and vectors used to compute tractions
        mfem::DenseMatrix F_, dNdeta_, dNdx_;
        mfem::Vector N_, T_vec_, tmp1_;

        // PK1 traction coefficient (not owned)
        mfem::VectorCoefficient* T_;

    public:
        
        PK1TractionIntegrator(mfem::VectorCoefficient& T) : T_(&T) {}

        void AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                    const mfem::Vector& elfun, mfem::Vector& elvect) override;

        void AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                    const mfem::Vector& elfun, mfem::DenseMatrix& elmat) override;

};

#endif
