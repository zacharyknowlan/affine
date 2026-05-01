#ifndef TL_STRESS_STRAIN
#define TL_STRESS_STRAIN

#include "mfem.hpp"

// Calculate the element-wise Green-Lagrange strain 
void CalcGreenLagrangeStrain(const mfem::GridFunction& u, mfem::GridFunction& E);

// Compute the element-wise Cauchy stress for the hyper elastic model
void CalcHyperElasticCauchyStress(const mfem::GridFunction& u, const mfem::GridFunction& E, 
                                        mfem::Coefficient& mu, mfem::Coefficient& lambda, 
                                        mfem::GridFunction& sigma);

// Compute the element-wise Cauchy stress for the Fung exponential model
void CalcFungCauchyStress(const mfem::GridFunction& u, const mfem::GridFunction& E, 
                            mfem::Coefficient& a, mfem::Coefficient& A1, mfem::Coefficient& A2, 
                            mfem::Coefficient& A3, mfem::Coefficient& A4, mfem::Coefficient& A5, 
                            mfem::Coefficient& A6, mfem::GridFunction& sigma);

// Compute the element-wise PK1 stress for the Fung exponential model
void CalcFungPK1Stress(const mfem::GridFunction& u, const mfem::GridFunction& E, 
                            mfem::Coefficient& a, mfem::Coefficient& A1, mfem::Coefficient& A2, 
                            mfem::Coefficient& A3, mfem::Coefficient& A4, mfem::Coefficient& A5, 
                            mfem::Coefficient& A6, mfem::GridFunction& P);

// Compute the element-wise Von-Mises stress
void CalcVonMisesStress(mfem::GridFunction& sigma, mfem::GridFunction& VMStress);

#endif 
