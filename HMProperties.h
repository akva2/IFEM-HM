// $Id$
//==============================================================================
//!
//! \file HMProperties.h
//!
//! \date Mar 24 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for properties for the heat-mass transfer problem.
//!
//==============================================================================

#ifndef _HM_PROPERTIES_H
#define _HM_PROPERTIES_H

#include "MatVec.h"
#include "Vec3.h"


namespace tinyxml2 { class XMLElement; }


/*!
  \brief Class representing properties for the heat-mass transfer problem.
*/

class HMProperties
{
public:
  //! \brief Empty constructor.
  HMProperties() = default;

  //! \brief Returns the elastic modulus as a function of temperature.
  //! \param T Temperature
  double elasticModulus(double T) const;

  //! \brief Parse settings from an xml element.
  //! \param elem The XML element to parse
  void parse(const tinyxml2::XMLElement* elem);
  //! \brief Print settings to terminal.
  void printLog();

  //! \brief Obtain the water velocity.
  //! \param T Temperature
  //! \param dC Concentration gradient
  //! \param dT Temperature gradient
  Vec3 waterVelocity(double T, const Vector& dC, const Vector& dT) const;

  //! \brief Returns water holding capacity as a function of temperature.
  //! \param T Temperature
  double waterCapacity(double T) const;

  //! \brief Returns the moisture diffusivity.
  double moistureDiffusivity() const { return diffusion; }
  //! \brief Returns the density of the meat.
  double meatDensity() const;
  //! \brief Returns the specific heat capacity of the meat.
  double meatHeatCapacity() const;
  //! \brief Returns the thermal conductivity of the meat.
  double meatThermalConductivity() const { return meat.conductivity; }
  //! \brief Returns the density of water.
  double waterDensity() const { return water.density; }
  //! \brief Returns the specific heat capacity of the water.
  double waterHeatCapacity() const { return water.capacity; }
  //! \brief Returns the latent vaporization heat of water.
  double vaporizationHeat() const { return water.vaporization_heat; }
  //! \brief Returns the fraction of energy that dissipates with evaporating water at the surface. 
  double dissipationCoefficient() const { return diss_coeff; }
  //! \brief Returns the oven temperature.
  double ovenTemperature() const { return T_oven; }
  //! \brief Returns the heat transfer coefficient.
  double heatTransfer() const { return heat_transfer; }

protected:
  //! \brief Returns the derivative of water holding capacity as a function of temperature.
  //! \param T Temperature
  double diffCeq(double T) const;
  //! \brief Returns the viscosity of water as a function of temperature.
  //! \param T Temperature
  double waterViscosity(double T) const;

  //! \brief Struct holding generic component parameters.
  struct Component {
    double density; //!< Density of component
    double composition; //!< Composition fraction for component
  };

  //! \brief Struct holding parameters for the meat.
  struct Meat : public Component {
    double conductivity; //!< Thermal conductivity of the meat
    double permeability; //!< Permeability of the meat
  };

  //! \brief Struct holding parameters for the water.
  struct Water : public Component {
    double capacity; //!< Heat capacity of water
    double vaporization_heat; //!< Latent vaporization heat of water
  };

  Component fat; //!< Properties of fat
  Component protein; //!< Properties of protein
  Component carbohydrate; //!< Properties of carbohydrates
  Water water; //!< Properties of water
  Meat meat; //!< Properties of meat
  double diffusion = 1.0; //!< Moisture diffusion coefficient
  double diss_coeff; //!< Dissipation coefficient
  double heat_transfer; //!< Heat transfer coefficient
  double vaporization_heat; //!< Latent vaporization heat of water
  double T_oven; //!< Oven temperature

  static constexpr double E0 = 12e3; //!< Minimum elastic modulus for meat (raw meat)
  static constexpr double Em = 83e3; //!< Maximum elastic modulus for meat (80C)
  static constexpr double En = 0.3; //!< Curve fitted for meat
  static constexpr double Ed = 60.0; //!< Curve fitted for meat
  static constexpr double Tsigma = 52; //!< Center of curve for excess water concentration
  static constexpr double a1 = 0.745; //!< Parameter in expression for excess water concentration
  static constexpr double a2 = 0.345; //!< Parameter in expression for excess water concentration
  static constexpr double a3 = 30; //!< Parameter in expression for excess water concentration
  static constexpr double a4 = 0.25; //!< Parameter in expression for excess water concentration
};

#endif
