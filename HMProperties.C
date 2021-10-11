// $Id$
//==============================================================================
//!
//! \file HMProperties.C
//!
//! \date Mar 24 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for properties for the heat-mass transferproblem.
//!
//==============================================================================

#include "HMProperties.h"
#include "IFEM.h"
#include "LogStream.h"
#include "Utilities.h"

#include "tinyxml.h"

#include <cmath>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>


void HMProperties::parse (const TiXmlElement* elem)
{
  auto&& parseComp = [elem](const char* name, Component& out)
  {
    const TiXmlElement* compElem = elem->FirstChildElement(name);
    if (compElem) {
      utl::getAttribute(compElem, "density", out.density);
      utl::getAttribute(compElem, "composition", out.composition);
    }
  };
  parseComp("fat", fat);
  parseComp("protein", protein);
  parseComp("carbohydrate", carbohydrate);
  parseComp("water", water);

  const TiXmlElement* waterComp = elem->FirstChildElement("water");
  if (waterComp) {
    utl::getAttribute(waterComp, "capacity", water.capacity);
    utl::getAttribute(waterComp, "vaporization_heat", water.vaporization_heat);
  }
  const TiXmlElement* meatComp = elem->FirstChildElement("meat");
  if (meatComp) {
    utl::getAttribute(meatComp, "conductivity", meat.conductivity);
    utl::getAttribute(meatComp, "permeability", meat.permeability);
  }

  utl::getAttribute(elem, "heat_transfer_coef", heat_transfer);
  utl::getAttribute(elem, "oven_temperature", T_oven);
  utl::getAttribute(elem, "diffusion", diffusion);
  utl::getAttribute(elem, "dissipation_coeff", diss_coeff);
}


void HMProperties::printLog ()
{
  IFEM::cout << "Heat-mass transfer problem properties:" << std::endl;
  IFEM::cout << "\t Heat transfer coefficient, h = " << heat_transfer << std::endl;
  IFEM::cout << "\t Oven temperature, To = " << T_oven << std::endl;
  IFEM::cout << "\t Diffusion coefficient, D = " << diffusion << std::endl;
  IFEM::cout << "\t Dissipation coefficient, f = " << diss_coeff << std::endl;
  auto&& printComp = [](const std::string& name, const Component& out)
  {
    IFEM::cout << "\t " << name <<" properties:" << std::endl;
    IFEM::cout << "\t\t Density, rho = " << out.density << std::endl;
    IFEM::cout << "\t\t Composition, y = " << out.composition << std::endl;
  };
  printComp("Fat", fat);
  printComp("Protein", protein);
  printComp("Carbohydrate", carbohydrate);
  printComp("Water", water);
  IFEM::cout << "\t\t Specific heat capacity, c = " << water.capacity << std::endl;
  IFEM::cout << "\t\t Latent heat of vaporization, H = " << water.vaporization_heat << std::endl;
  IFEM::cout << "\tProperties of meat:" << std::endl;
  IFEM::cout << "\t\t Density, rho = " << this->meatDensity() << std::endl;
  IFEM::cout << "\t\t Specific heat capacity, c = " << this->meatHeatCapacity() << std::endl;
  IFEM::cout << "\t\t Thermal conductivity, k = " << meat.conductivity << std::endl;
  IFEM::cout << "\t\t Permeability, K = " << meat.permeability << std::endl;
}


double HMProperties::waterCapacity (double T) const
{
  return a1 - a2 / (1.0 + a3*exp(-a4*(T-Tsigma))); // (7)
}


double HMProperties::elasticModulus (double T) const
{
  return E0 + Em / (1.0 + exp(-En*(T-Ed))); // (8)
}


double HMProperties::diffCeq (double T) const
{
  return -a2*a3*a4*exp(a4*(T-Tsigma)) / pow(a3 + exp(a4*(T-Tsigma)), 2.0); // derivative of (7)
}


Vec3 HMProperties::waterVelocity (double T,
                                  const Vector& dC, const Vector& dT) const
{
  double scale = -meat.permeability * elasticModulus(T) / this->waterViscosity(T);
  double dCeq = diffCeq(T);

  Vec3 result;
  for (size_t i = 0; i < 3; ++i)
      result[i] = scale * (dC[i] - dCeq*dT[i]);

  return result;
}


double HMProperties::meatDensity () const
{
  double denom = fat.composition / fat.density + 
                 protein.composition / protein.density +
                 carbohydrate.composition / carbohydrate.density +
                 water.composition / water.density;

  return 1.0 / denom;
}


double HMProperties::meatHeatCapacity () const
{
  return (1.6*carbohydrate.composition + 2.0*protein.composition +
          2.0*fat.composition + 4.2*water.composition) * 1e3;
}


double HMProperties::waterViscosity (double T) const
{
  return exp(-0.0072*T - 2.8658);
}
