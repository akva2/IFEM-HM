// $Id$
//==============================================================================
//!
//! \file MassTransfer.C
//!
//! \date Mar 25 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for mass transfer problems.
//!
//==============================================================================

#include "MassTransfer.h"

#include "ElmMats.h"
#include "FiniteElement.h"
#include "LocalIntegral.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"

#include <ext/alloc_traits.h>
#include <iostream>
#include <memory>


MassTransfer::MassTransfer (unsigned short int n,
                            TimeIntegration::Method method) :
  IntegrandBase(n),
  bdf(TimeIntegration::Order(method))
{
  primsol.resize(TimeIntegration::Steps(method) + 1);
  registerVector("temperature1",&temperature);
}


bool MassTransfer::initElement (const std::vector<int>& MNPC,
                                LocalIntegral& A)
{
  if (!this->IntegrandBase::initElement(MNPC, A))
      return false;

  A.vec.resize(A.vec.size()+1);
  int ierr;
  if ((ierr = utl::gather(MNPC,1,temperature,A.vec.back())))
    std::cerr <<" *** " << __PRETTY_FUNCTION__ << ":  Detected "
              << ierr <<" node numbers out of range."<< std::endl;

  return ierr == 0;
}


bool MassTransfer::evalInt (LocalIntegral& elmInt,
                            const FiniteElement& fe,
                            const TimeDomain& time,
                            const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  double Cbdf = 0;
  for (int i = 1; i <= bdf.getOrder(); ++i)
    Cbdf += -fe.N.dot(elMat.vec[i]) * bdf[i] / time.dt;

  double T = elMat.vec.back().dot(fe.N);
  Vector dC;
  fe.dNdX.multiply(elMat.vec[1], dC, true);
  Vector dT;
  fe.dNdX.multiply(elMat.vec.back(), dT, true);

  WeakOps::Mass(elMat.A[0], fe, bdf[0] / time.dt);
  WeakOps::Laplacian(elMat.A[0], fe, props.moistureDiffusivity());

  Vec3 Uw = props.waterVelocity(T, dC, dT);
  WeakOps::Advection(elMat.A[0], fe, Uw, 1.0, WeakOperators::CONSERVATIVE);

  WeakOps::Source(elMat.b[0], fe, Cbdf);

  return true;
}


std::string MassTransfer::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "C";

  return prefix + std::string(" C");
}


MassTransfer::Robin::Robin(unsigned short int n, const MassTransfer& itg) :
  IntegrandBase(n),
  integrand(itg)
{
}


bool MassTransfer::Robin::initElementBou(const std::vector<int>& MNPC,
                                                 LocalIntegral& elmInt)
{
  return const_cast<MassTransfer&>(integrand).initElement(MNPC, elmInt);
}


bool MassTransfer::Robin::evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                                          const Vec3& X, const Vec3& normal) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  const HMProperties& props = integrand.getProps();

  double C = elMat.vec[1].dot(fe.N);
  double T = elMat.vec.back().dot(fe.N);
  Vector dT;
  fe.dNdX.multiply(elMat.vec.back(), dT, true);
  Vector dC;
  fe.dNdX.multiply(elMat.vec[1], dC, true);

  Vec3 uw = props.waterVelocity(T, dC, dT);

  double delT = props.ovenTemperature() - T;
  double denom = props.vaporizationHeat() * props.waterDensity();
  double val = props.dissipationCoefficient() * props.heatTransfer() * delT / denom;

  WeakOps::Mass(elMat.A[0], fe, normal*uw);
  WeakOps::Source(elMat.b[0], fe, val * (C - props.waterCapacity(T)));

  return true;
}
