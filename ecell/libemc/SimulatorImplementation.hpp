//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-CELL Simulation Environment package
//
//                Copyright (C) 1996-2000 Keio University
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-CELL is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// E-CELL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with E-CELL -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Kouichi Takahashi <shafi@e-cell.org> at
// E-CELL Project, Lab. for Bioinformatics, Keio University.
//


#ifndef ___SIMULATOR_IMPLEMENTATION_H___
#define ___SIMULATOR_IMPLEMENTATION_H___

#include "libecs.hpp"

/**
   Pure virtual base class (interface definition) of simulator
   implementation.
*/

class SimulatorImplementation
{

public:

  SimulatorImplementation() {}
  virtual ~SimulatorImplementation() {}

  virtual RootSystemRef getRootSystem() = 0;

  virtual void createEntity( StringCref classname, 
			      FQPICref fqpi, 
			      StringCref name ) = 0;

  virtual void setProperty( FQPICref fqpi, 
			    MessageCref message ) = 0;

  virtual const Message getProperty( FQPICref fqpi, 
			      StringCref propertyName ) = 0;
  virtual void step() = 0;

  virtual void initialize() = 0;

};   //end of class Simulator

#endif   /* ___SIMULATOR_IMPLEMENTATION_H___ */













