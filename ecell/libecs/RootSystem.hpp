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

#ifndef ___ROOTSYSTEM_H___
#define ___ROOTSYSTEM_H___
#include "System.hpp"
#include "SubstanceMaker.hpp"
#include "ReactorMaker.hpp"
#include "SystemMaker.hpp"
#include "AccumulatorMaker.hpp"


class RootSystem : public System
{

public:

  class MalformedSystemName : public NotFound
    {
    public:
      MalformedSystemName( StringCref method, StringCref message ) 
	: NotFound( method, message ) {}
      const String what() const { return "Malformed system name."; }
    };

public:

  RootSystem();
  ~RootSystem();

  int check();

  SystemPtr getSystem( SystemPathCref systempath )
    throw( NotFound, MalformedSystemName );
  Entity getEntity( FQPICref fqpi ) 
    throw( InvalidPrimitiveType, NotFound );

  virtual void initialize();

  StepperLeaderRef    getStepperLeader()    { return theStepperLeader; }

  ReactorMakerRef     getReactorMaker()     { return theReactorMaker; }
  SubstanceMakerRef   getSubstanceMaker()   { return theSubstanceMaker; }
  SystemMakerRef      getSystemMaker()      { return theSystemMaker; }
  StepperMakerRef     getStepperMaker()     { return theStepperMaker; }
  AccumulatorMakerRef getAccumulatorMaker() { return theAccumulatorMaker; }

  virtual const char* const className() const { return "RootSystem"; }

private:

  void install();

private:

  StepperLeader    theStepperLeader;

  ReactorMaker     theReactorMaker;
  SubstanceMaker   theSubstanceMaker;
  SystemMaker      theSystemMaker;
  StepperMaker     theStepperMaker;
  AccumulatorMaker theAccumulatorMaker;

};

#endif /* ___ROOTSYSTEM_H___ */


/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/
