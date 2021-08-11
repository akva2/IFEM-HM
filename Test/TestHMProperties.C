//==============================================================================
//!
//! \file TestHMProperties.C
//!
//! \date Aug 11 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for driver for heat-mass transfer properties.
//!
//==============================================================================

#include "HMProperties.h"

#include "gtest/gtest.h"
#include "tinyxml.h"

TEST(TestHMProperties, Parse)
{
  HMProperties hm;
  TiXmlDocument doc;
  doc.Parse(R"(<hmproperties heat_transfer_coef="33.4" oven_temperature="175"
              diffusion="4e-10" dissipation_coeff="0.88">
              <fat density="920" composition="0.03"/>
              <protein density="1320" composition="0.2"/>
              <carbohydrate density="1600" composition="0.02"/>
              <water density="998" composition="0.75" capacity="4170" vaporization_heat="2.3e6"/>
              <meat conductivity="0.4" permeability="1e-17"/>
              </hmproperties>)");
  hm.parse(doc.RootElement());

  EXPECT_FLOAT_EQ(hm.meatDensity(), 1054.7112);
  EXPECT_FLOAT_EQ(hm.meatHeatCapacity(), 3642.0);
  EXPECT_FLOAT_EQ(hm.meatThermalConductivity(), 0.4);
}

