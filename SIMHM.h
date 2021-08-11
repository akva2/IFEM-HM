// $Id$
//==============================================================================
//!
//! \file SIMHM.h
//!
//! \date Mar 25 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for heat-mass transfer simulators.
//!
//==============================================================================

#ifndef _SIM_HM_H_
#define _SIM_HM_H_

#include "TimeStep.h"
#include "SIMCoupledSI.h"
#include "SIMsolution.h"


/*!
  \brief Driver class for segregated heat-mass transfer simulators.
*/

template<class HeatSolver, class MassSolver,
         template<class T1, class T2> class Coupling>
class SIMHMBase : public Coupling<HeatSolver,MassSolver>
{
public:
  //! \brief Enum announcing the dimensionality (used for template writing).
  enum { dimension = HeatSolver::dimension };

  //! \brief The constructor initializes the references to the two solvers.
  SIMHMBase(HeatSolver& s1, MassSolver& s2)
    : Coupling<HeatSolver,MassSolver>(s1,s2)
  {
  }

  //! \brief Empty destructor. We don't own our subsolvers.
  virtual ~SIMHMBase() {}

  //! \brief Initializes and sets up field dependencies.
  void setupDependencies() override
  {
    this->S1.registerDependency(&this->S2, "concentration1");
    this->S2.registerDependency(&this->S1, "temperature1");
  }
};


template<class HeatSolver, class MassSolver>
using SIMHM = SIMHMBase<HeatSolver,MassSolver,SIMCoupled>; //!< Convenience typedef


/*!
  \brief Driver class for semi-implicit heat-mass transfer simulators.
*/

template<class HeatSolver, class MassSolver>
class SIMHMSI : public SIMHMBase<HeatSolver, MassSolver, SIMCoupledSI>
{
public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMHMSI(HeatSolver& s1, MassSolver& s2)
    : SIMHMBase<HeatSolver, MassSolver, SIMCoupledSI>(s1, s2)
  {
  }

  //! \brief Check convergence of the total system.
  SIM::ConvStatus checkConvergence(const TimeStep& tp,
                                   SIM::ConvStatus status1,
                                   SIM::ConvStatus status2) override
  {
    if (status1 == SIM::FAILURE || status2 == SIM::FAILURE)
      return SIM::FAILURE;
    else if (status1 == SIM::DIVERGED || status2 == SIM::DIVERGED)
      return SIM::DIVERGED;

    double rConv = 1.0;
    if (!prevTemp.empty()) {
      Vector dTemp = static_cast<const SIMsolution&>(this->S1).getSolution(0);
      dTemp -= prevTemp;
      Vector dMass = static_cast<const SIMsolution&>(this->S2).getSolution(0);
      dMass -= prevMass;
      rConv = sqrt(pow(dTemp.norm2(), 2.0) + pow(dMass.norm2(), 2.0));
    }
    prevTemp = static_cast<const SIMsolution&>(this->S1).getSolution(0);
    prevMass = static_cast<const SIMsolution&>(this->S2).getSolution(0);

    int maxCycle = this->S1.getMaxit();
    double cycleTol = this->S1.getSubItTol();
    IFEM::cout <<"  cycle "<< tp.iter <<": Res = "<< rConv << std::endl;
    if (rConv < cycleTol) {
      prevTemp.clear();
      prevMass.clear();
      return SIM::CONVERGED;
    } else if (tp.iter < maxCycle)
      return SIM::OK;

    std::cerr <<"SIMBeefSI::checkConvergence: Did not converge in "
              << maxCycle <<" staggering cycles, bailing.."<< std::endl;
    return SIM::DIVERGED;
  }

protected:
  Vector prevTemp; //!< Previous temperature solution
  Vector prevMass; //!< Previous mass solution
};

#endif
