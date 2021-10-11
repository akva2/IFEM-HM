// $Id$
//==============================================================================
//!
//! \file main_HeatTransfer.C
//!
//! \date Mar 13 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Main program for the isogeometric heat transfer solver.
//!
//==============================================================================

#include "HeatTransfer.h"
#include "HMArgs.h"
#include "SIMHeatTransfer.h"

#include "IFEM.h"
#include "LogStream.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMadmin.h"
#include "SIMoptions.h"
#include "Profiler.h"
#include "Property.h"
#include "SIMSolver.h"
#include "TimeStep.h"

#include <iostream>
#include <map>
#include <vector>


/*!
  \brief Runs the transient beef problem.
*/

template<class Dim>
int runSimulator (char* infile, const HMArgs& args)
{
  using HeatSim = SIMHeatTransfer<Dim, HeatTransfer>;

  utl::profiler->start("Model input");
  HeatTransfer heat(Dim::dimension, false, args.timeMethod);
  HeatSim heatSim(heat);

  SIMSolver<HeatSim> solver(heatSim);

  // Read in model definitions
  if (!heatSim.read(infile))
    return 1;

  // Read in model definitions
  if (!solver.read(infile))
    return 1;

  if (!heatSim.preprocess())
    return 2;

  if (!heatSim.init(TimeStep()))
    return 2;

  heatSim.setInitialConditions();

  heatSim.opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  if (solver.restart(heatSim.opt.restartFile,heatSim.opt.restartStep) < 0)
    return 2;

  if (heatSim.opt.dumpHDF5(infile))
    solver.handleDataOutput(heatSim.opt.hdf5,heatSim.getProcessAdm(),
                            heatSim.opt.saveInc,heatSim.opt.restartInc);

  int res = solver.solveProblem(infile,"Solving Heat transfer problem");
  return res;
}


/*!
  \brief Main program for the isogeometric Heat transfer solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -2D : Use two-parametric simulation driver
  \arg -adap : Use adaptive simulation driver with LR-splines discretization
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  char* infile = nullptr;
  IFEM::Init(argc,argv,"Heat transfer solver");
  HMArgs args;
  for (int i = 1; i < argc; i++)
    if (argv[i] == infile || args.parseArg(argv[i]))
      ; // ignore the input file on second pass
    else if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!infile) {
      infile = argv[i];
      if (!args.readXML(infile,false))
        return 1;
      i = 0;
    }
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-2D] [-nGauss <n>] [-adap]\n"
              <<"       [-hdf5] [-vtf <format> [-nviz <nviz>]"
              <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]]\n";
    return 0;
  }

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout) << std::endl;
  utl::profiler->stop("Initialization");
  args.print();

  SIMadmin::msgLevel = 1; // Disable verbose output

  if (args.dim == 2)
    return runSimulator<SIM2D>(infile, args);
  else
    return runSimulator<SIM3D>(infile, args);
}
