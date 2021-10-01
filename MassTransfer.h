// $Id$
//==============================================================================
//!
//! \file MassTransfer.h
//!
//! \date Mar 13 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for mass transfer in Beef problem.
//!
//==============================================================================

#ifndef _MASS_TRANSFER_H
#define _MASS_TRANSFER_H

#include "BDF.h"
#include "EqualOrderOperators.h"
#include "IntegrandBase.h"
#include "TimeIntUtils.h"

#include "HMProperties.h"


/*!
  \brief Class representing the integrand of mass transfer in a beef problem.
*/

class MassTransfer : public IntegrandBase
{
public:
  using WeakOps = EqualOrderOperators::Weak; //!< Convenience alias

  //! \brief Class representing the weak Dirichlet integrand.
  class Robin : public IntegrandBase
  {
  public:
    //! \brief Default constructor.
    //! \param[in] n Number of spatial dimensions
    //! \param[in] itg Main integrand instance
    Robin(unsigned short int n, const MassTransfer& itg);
    //! \brief Empty destructor.
    virtual ~Robin() {}

    //! \brief Returns that this integrand has no interior contributions.
    bool hasInteriorTerms() const override { return false; }

    using IntegrandBase::getLocalIntegral;
    //! \brief Returns a local integral contribution object for given element.
    //! \param[in] nen Number of nodes on element
    //! \param[in] iEl Element index
    LocalIntegral* getLocalIntegral(size_t nen, size_t iEl, bool) const override
    {
      return integrand.getLocalIntegral(nen, iEl, false);
    }

    using IntegrandBase::initElementBou;
    //! \brief Initializes current element for boundary integration.
    //! \param[in] MNPC Matrix of nodal point correspondance for current element
    //! \param elmInt Local integral for element
    bool initElementBou(const std::vector<int>& MNPC,
                        LocalIntegral& elmInt) override;

    using IntegrandBase::evalBou;
    //! \brief Evaluates the integrand at a boundary point.
    //! \param elmInt The local integral object to receive the contributions
    //! \param[in] fe Finite element data of current integration point
    //! \param[in] X Cartesian coordinates of current integration point
    //! \param[in] normal Boundary normal vector at current integration point
    bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                 const Vec3& X, const Vec3& normal) const override;

  protected:
    const MassTransfer& integrand; //!< Main integrand instance
  };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] method The time integration method to use
  explicit MassTransfer(unsigned short int,
                        TimeIntegration::Method method = TimeIntegration::BE);

  //! \brief Empty destructor.
  virtual ~MassTransfer() {}

  using IntegrandBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param A The local integral to initialize
  bool initElement(const std::vector<int>& MNPC, LocalIntegral& A) override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

   //! \brief Returns a const ref to the problem properties.
   const HMProperties& getProps() const { return props; }
   //! \brief Returns a reference to the problem properties.
   HMProperties& getProps() { return props; }

   //! \brief Advance time stepping
   void advanceStep() { bdf.advanceStep(); }

   //! \brief Returns the name of the primary solution field.
   //! \param[in] prefix Name prefix
   std::string getField1Name (size_t, const char* prefix) const override;

protected:
  HMProperties props; //!< Material properties
  Vector temperature; //!< Externally provided temperature
  TimeIntegration::BDF bdf; //!< BDF time stepping helper
};

#endif
