// $Id$
//==============================================================================
//!
//! \file HeatTransfer.C
//!
//! \date Mar 25 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for heat transfer problems.
//!
//==============================================================================

#include "HeatTransfer.h"
#include "ElmMats.h"
#include "FiniteElement.h"
#include "LocalIntegral.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3Oper.h"

#include <iostream>
#include <memory>


HeatTransfer::HeatTransfer (unsigned short int n,
                            bool mass,
                            TimeIntegration::Method method) :
  IntegrandBase(n),
  withMass(mass),
  bdf(TimeIntegration::Order(method))
{
  primsol.resize(TimeIntegration::Steps(method) + 1);
  registerVector("concentration1",&concentration);
}


bool HeatTransfer::initElement (const std::vector<int>& MNPC,
                                LocalIntegral& A)
{
  if (!this->IntegrandBase::initElement(MNPC, A))
      return false;

  int ierr = 0;
  if (withMass) {
    A.vec.resize(A.vec.size()+1);
    if ((ierr = utl::gather(MNPC,1,concentration,A.vec.back())))
      std::cerr <<" *** " << __PRETTY_FUNCTION__ << ":  Detected "
                << ierr <<" node numbers out of range."<< std::endl;
  }

  return ierr == 0;
}


bool HeatTransfer::evalInt (LocalIntegral& elmInt,
                            const FiniteElement& fe,
                            const TimeDomain& time,
                            const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  double T = elMat.vec[1].dot(fe.N);
  double Tbdf = -T * bdf[1] / time.dt;
  for (int i = 2; i <= bdf.getOrder(); ++i)
    Tbdf -= fe.N.dot(elMat.vec[i]) * bdf[i] / time.dt;

  WeakOps::Mass(elMat.A[0], fe, props.meatDensity() * props.meatHeatCapacity() * bdf[0] / time.dt);
  WeakOps::Laplacian(elMat.A[0], fe, props.meatThermalConductivity());

  if (withMass) {
    Vector dT;
    fe.dNdX.multiply(elMat.vec[1], dT, true);
    Vector dC;
    fe.dNdX.multiply(elMat.vec.back(), dC, true);
    Vec3 waterVelocity = props.waterVelocity(T, dC, dT);
    WeakOps::Advection(elMat.A[0], fe, waterVelocity,
                       props.waterDensity() * props.waterHeatCapacity(),
                       WeakOperators::CONSERVATIVE);
  }

  WeakOps::Source(elMat.b[0], fe, props.meatDensity() * props.meatHeatCapacity() * Tbdf);

  return true;
}


bool HeatTransfer::evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                           const Vec3& X, const Vec3& normal) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  double T = elMat.vec[1].dot(fe.N);
  WeakOps::Source(elMat.b[0], fe, props.heatTransfer() * (props.ovenTemperature() - T));

  return true;
}


std::string HeatTransfer::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "T";

  return prefix + std::string(" T");
}


HeatTransfer::Robin::Robin(unsigned short int n, const HeatTransfer& itg) :
  IntegrandBase(n),
  integrand(itg)
{
}


bool HeatTransfer::Robin::initElementBou(const std::vector<int>& MNPC,
                                                 LocalIntegral& elmInt)
{
  return const_cast<HeatTransfer&>(integrand).initElement(MNPC, elmInt);
}


bool HeatTransfer::Robin::evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                                          const Vec3& X, const Vec3& normal) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  const HMProperties& props = integrand.getProps();

  if (integrand.enableMass()) {
    double T = elMat.vec[1].dot(fe.N);
    Vector dT;
    fe.dNdX.multiply(elMat.vec[1], dT, true);
    Vector dC;
    fe.dNdX.multiply(elMat.vec.back(), dC, true);
    Vec3 uw = props.waterVelocity(T, dC, dT);
    double val = (1.0 - props.dissipationCoefficient()) * props.heatTransfer();
    WeakOps::Mass(elMat.A[0], fe, (uw * normal) * props.waterHeatCapacity() * props.waterDensity());
    WeakOps::Source(elMat.b[0], fe, val * (props.ovenTemperature() - T));
  } else {
    std::cerr << "No Robin conditions without mass contributions, use Neumann conditions." << std::endl;
    return false;
  }

  return true;
}
