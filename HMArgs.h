// $Id$
//==============================================================================
//!
//! \file HMArgs.h
//!
//! \date Mar 25 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the heat-mass transfer application.
//!
//==============================================================================

#ifndef _HM_ARGS_H
#define _HM_ARGS_H

#include "SIMargsBase.h"
#include "Integrand.h"
#include "TimeIntUtils.h"


/*!
  \brief Class holding application parameters.
*/

class HMArgs : public SIMargsBase
{
public:
  TimeIntegration::Method timeMethod = TimeIntegration::BE; //!< Time integration method
  bool semiImplicit = false; //!< To use sub-iterations

  //! \brief Default constructor.
  HMArgs() : SIMargsBase("heatmasstransfer") {}

  //! \brief Parses a command-line argument.
  bool parseArg(const char* argv) override;

  //! \brief Print settings to terminal.
  void print();

protected:
  //! \brief Parse an element from the input file.
  bool parse(const TiXmlElement* elem) override;
};

#endif
