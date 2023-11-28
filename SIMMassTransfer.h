// $Id$
//==============================================================================
//!
//! \file SIMMassTransfer.h
//!
//! \date Mar 13 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for the mass transfer in the Beef problem.
//!
//==============================================================================

#ifndef _SIM_MASS_TRANSFER_H
#define _SIM_MASS_TRANSFER_H

#include "DataExporter.h"
#include "Property.h"
#include "Profiler.h"
#include "SIMenums.h"
#include "SIMMultiPatchModelGen.h"
#include "SIMsolution.h"
#include "TimeStep.h"
#include "tinyxml2.h"
#include "IFEM.h"


/*!
  \brief Driver class for the mass transfer in the Beef problem simulator.
*/

template<class Dim, class Integrand>
class SIMMassTransfer : public SIMMultiPatchModelGen<Dim>, public SIMsolution
{
public:
  //! \brief Default constructor.
  //! \param[in] m Integrand for mass transfer problem
  explicit SIMMassTransfer(Integrand& m) :
    SIMMultiPatchModelGen<Dim>(1),
    mass(m),
    robinBC(Dim::dimension, mass)
  {
    Dim::myHeading = "Mass transfer solver";
    Dim::myProblem = &mass;
  }

  virtual ~SIMMassTransfer()
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "MassTransfer"; }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override
  {
    if (!strcasecmp(elem->Value(), "hmproperties")) {
      mass.getProps().parse(elem);    
      return true;
    }

    if (strcasecmp(elem->Value(),"masstransfer"))
      return this->Dim::parse(elem);

    const tinyxml2::XMLElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      this->Dim::parse(child);

    return true;
  }

  //! \brief Initialize simulator.
  bool init(const TimeStep&)
  {
    // Initialize temperature solution vectors
    size_t n, nSols = this->getNoSolutions();
    this->initSolution(this->getNoDOFs(),nSols);
    std::string str = "concentration1";
    for (n = 0; n < nSols; n++, str[13]++)
      this->registerField(str,solution[n]);

    if (!this->initSystem(Dim::opt.solver))
      return false;

    this->setQuadratureRule(Dim::opt.nGauss[0],true);

    return true;
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    PROFILE1("SIMMassTransfer::solveStep");

    this->setMode(SIM::DYNAMIC);
    Vector dummy;
    this->updateDirichlet(tp.time.t,&dummy);

    if (!this->assembleSystem(tp.time,solution))
      return false;

    if (!this->solveSystem(solution.front(),0,"concentration "))
      return false;

    double dMax;
    size_t iMax;
    double normL2 = this->solutionNorms(solution.front(),&dMax,&iMax,1);
    IFEM::cout <<"\n  Concentration summary:  L2-norm        : "<< normL2
               <<"\n                       Max concentration : "<< dMax << " node " << iMax
               << std::endl;

    return true;
  }

  //! \brief No post-processing.
  bool postSolve(TimeStep&) { return true; }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep&)
  {
    this->pushSolution(); // Update solution vectors between time steps
    //mass.advanceStep();
    return true;
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMMassTransfer::saveStep");

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    if (!this->writeGlvS(solution.front(),iDump,nBlock,
                         tp.time.t,"concentration",99))
      return false;

    return this->writeGlvStep(iDump,tp.time.t);
  }

  //! \brief Register fields for simulation result export.
  void registerFields(DataExporter& exporter)
  {
    int results = DataExporter::PRIMARY;
    exporter.registerField("T","temperature",DataExporter::SIM,results);
    exporter.setFieldValue("T", this, &solution.front());
  }

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(SerializeMap& data) const override
  {
    return this->saveSolution(data,this->getName());
  }

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const SerializeMap& data) override
  {
    if (!this->restoreSolution(data,this->getName()))
      return false;

    mass.advanceStep();
    return true;
  }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //! \details Since this solver is linear, this is just a normal solve.
  SIM::ConvStatus solveIteration(TimeStep& tp)
  {
    return this->solveStep(tp) ? SIM::CONVERGED : SIM::FAILURE;
  }

  //! \brief Preprocessing performed before the FEM model generation.
  //! \details This method is reimplemented to couple the weak Dirichlet
  //! integrand to the Robin property codes.
  void preprocessA() override
  {
    Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

    for (const Property& p : Dim::myProps)
      if (p.pcode == Property::ROBIN)
        if (Dim::myInts.find(p.pindx) == Dim::myInts.end())
          Dim::myInts.insert(std::make_pair(p.pindx,&robinBC));
  }

protected:
  Integrand& mass; //!< Main integrand
  typename Integrand::Robin robinBC; //!< Integrand for Robin boundary condition
};

#endif
