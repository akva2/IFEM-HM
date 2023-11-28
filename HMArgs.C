// $Id$
//==============================================================================
//!
//! \file HMArgs.C
//!
//! \date Mar 25 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the heat-mass transfer application.
//!
//==============================================================================

#include "HMArgs.h"
#include "IFEM.h"
#include "LogStream.h"
#include "Utilities.h"
#include "tinyxml2.h"

#include <iostream>
#include <string>
#include <strings.h>


bool HMArgs::parseArg (const char* argv)
{
  TimeIntegration::Method tmp;
  if (argv[0] != '-')
    return false;
  else if ((tmp = TimeIntegration::get(argv+1)) > TimeIntegration::NONE) {
    if (tmp!= TimeIntegration::BE && tmp != TimeIntegration::BDF2)
      tmp = TimeIntegration::BE;
    timeMethod = tmp;
  } else if (!strcasecmp(argv,"-semiimplicit"))
    semiImplicit = true;
  else
    return this->SIMargsBase::parseArg(argv);

  return true;
}


bool HMArgs::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"timestepping")) {
    std::string type;
    if (utl::getAttribute(elem,"type",type))
      timeMethod = TimeIntegration::get(type);
    if (timeMethod != TimeIntegration::BE && timeMethod != TimeIntegration::BDF2)
      timeMethod = TimeIntegration::BE;
  }

  return this->SIMargsBase::parse(elem);
}


void HMArgs::print ()
{
  IFEM::cout << "Time integration: ";
  if (timeMethod == TimeIntegration::BE)
    IFEM::cout << "Backward Euler";
  else
    IFEM::cout << "BDF2";
  IFEM::cout << std::endl;
  if (semiImplicit)
    IFEM::cout << "Sub-iterations are enabled." << std::endl;
}
