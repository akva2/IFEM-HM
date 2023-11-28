// $Id$
//==============================================================================
//!
//! \file SIMHeatTransfer.h
//!
//! \date Mar 13 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for the heat transfer in the Beef problem.
//!
//==============================================================================

#ifndef _SIM_HEAT_TRANSFER_H
#define _SIM_HEAT_TRANSFER_H

#include "DataExporter.h"
#include "IFEM.h"
#include "Profiler.h"
#include "Property.h"
#include "SIMenums.h"
#include "SIMMultiPatchModelGen.h"
#include "SIMsolution.h"
#include "TimeStep.h"

#include "tinyxml2.h"


/*!
  \brief Driver class for the heat transfer in the Beef problem simulator.
*/

template<class Dim, class Integrand>
class SIMHeatTransfer : public SIMMultiPatchModelGen<Dim>, public SIMsolution
{
public:
  //! \brief Default constructor.
  //! \param[in] h Integrand for heat transfer problem
  explicit SIMHeatTransfer(Integrand& h) :
    SIMMultiPatchModelGen<Dim>(1),
    heat(h),
    robinBC(Dim::dimension, heat)
  {
    Dim::myHeading = "Heat transfer solver";
    Dim::myProblem = &heat;
  }

  virtual ~SIMHeatTransfer() 
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "HeatTransfer"; }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override
  {
    if (!strcasecmp(elem->Value(), "hmproperties")) {
      heat.getProps().parse(elem);    
      heat.getProps().printLog();
      return true;
    }

    if (strcasecmp(elem->Value(),"heattransfer"))
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
    std::string str = "temperature1";
    for (n = 0; n < nSols; n++, str[11]++)
      this->registerField(str,solution[n]);

    if (!this->initSystem(Dim::opt.solver))
      return false;

    this->setQuadratureRule(Dim::opt.nGauss[0],true);

    return true;
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    PROFILE1("SIMHeatTransfer::solveStep");

    this->setMode(SIM::DYNAMIC);
    if (Dim::msgLevel >= 0 && (!subIter || tp.iter == 0))
      IFEM::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    Vector dummy;
    this->updateDirichlet(tp.time.t,&dummy);

    if (!this->assembleSystem(tp.time,solution))
      return false;

    if (!this->solveSystem(solution.front(),0,"temperature "))
      return false;

    double dMax;
    size_t iMax;
    double normL2 = this->solutionNorms(solution.front(),&dMax,&iMax,1);
    IFEM::cout <<"\n  Temperature summary:  L2-norm        : "<< normL2
               <<"\n                       Max temperature : "<< dMax << " node " << iMax
               << std::endl;

    return true;
  }

  //! \brief No post-processing.
  bool postSolve(TimeStep&) { return true; }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep&)
  {
    this->pushSolution(); // Update solution vectors between time steps
    heat.advanceStep();
    return true;
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0) return true;

    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMHeatTransfer::saveStep");

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    if (!this->writeGlvS(solution.front(),iDump,nBlock,
                         tp.time.t,"temperature",89))
      return false;

    return heat.enableMass() ? true : this->writeGlvStep(iDump,tp.time.t);
  }

  //! \brief Register fields for simulation result export.
  void registerFields(DataExporter& exporter)
  {
    int results = DataExporter::PRIMARY;
    exporter.registerField("C","concentration",DataExporter::SIM,results);
    exporter.setFieldValue("C", this, &solution.front());
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

    heat.advanceStep();
    return true;
  }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //! \details Since this solver is linear, this is just a normal solve.
  SIM::ConvStatus solveIteration(TimeStep& tp)
  {
    subIter = true;
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
  bool subIter = false; //!< True if we are doing sub-iterations
  Integrand& heat; //!< Main integrand
  typename Integrand::Robin robinBC; //!< Integrand for Robin boundary conditions 
};

#endif
