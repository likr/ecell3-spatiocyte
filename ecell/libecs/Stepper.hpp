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

#ifndef ___STEPPER_H___
#define ___STEPPER_H___
#include "libecs.hpp"

typedef IntegratorPtr ( *IntegratorAllocator_ )( Substance& );
DECLARE_TYPE( IntegratorAllocator_, IntegratorAllocator );

typedef Stepper (* StepperAllocatorFunc )();


class Stepper
{

public:

  Stepper(); 
  virtual ~Stepper() {}

  void setOwner( SystemPtr owner ) { theOwner = owner; }
  SystemPtr getOwner() const { return theOwner; }

  virtual void initialize();

  virtual Float getDeltaT() = 0;
  virtual bool isMasterStepper() const { return false; }

  virtual void clear() = 0;
  virtual void react() = 0;
  virtual void transit() = 0;
  virtual void postern() = 0;

  virtual void distributeIntegrator( IntegratorAllocatorPtr );

  virtual const char* const className() const  { return "Stepper"; }

protected:

  SystemPtr theOwner;

};

class MasterStepper : public Stepper
{

public:

  typedef list<SlaveStepperPtr> SlaveStepperList;
  typedef SlaveStepperList::iterator SlaveStepperListIterator;

  MasterStepper();

  void setPace( int pace ) { thePace = pace; }
  void setdeltaT( Float dt ) { theDeltaT = dt; }

  virtual ~MasterStepper() {}

  virtual int getNumberOfSteps() const { return 1; }
  int getPace() const { return thePace; }

  virtual Float getDeltaT(); 
  virtual void initialize();

  virtual void distributeIntegrator( IntegratorAllocator );
  void registerSlaves( SystemPtr );

  virtual bool isMasterStepper() const { return true; }

  virtual const char* const className() const  { return "MasterStepper"; }

protected:

  int                 thePace;
  Float               theDeltaT;
  IntegratorAllocator theAllocator;
  SlaveStepperList    theSlavesList;

};


class StepperLeader
{

  typedef multimap<int,MasterStepperPtr> MasterStepperMap;

public:

  StepperLeader();
  virtual ~StepperLeader() {}

  void registerMasterStepper( MasterStepperPtr newstepper );

  static void setDefaultUpdateDepth( int d ) { DEFAULT_UPDATE_DEPTH = d; }
  static int getDefaultUpdateDepth()         { return DEFAULT_UPDATE_DEPTH; }

  void setUpdateDepth( int d ) { theUpdateDepth = d; }
  int getUpdateDepth()         { return theUpdateDepth; }

  virtual void initialize();

  Float getDeltaT();
  int getBaseClock() { return theBaseClock; }
  void step();
  virtual void clear();
  virtual void react();
  virtual void transit();
  virtual void postern();

  void update();
  virtual const char* const className() const  { return "StepperLeader"; }

protected:

  void setBaseClock( int baseclock );

protected:

  MasterStepperMap theStepperList;

private:

  int theUpdateDepth;
  int theBaseClock;

  static int DEFAULT_UPDATE_DEPTH;

};


#include "System.hpp"

class SlaveStepper : public Stepper
{

public:

  SlaveStepper() {}
  virtual ~SlaveStepper() {}

  
  // virtual and inlined.
  // this is for optimization in operations with SlaveStepperPtrs.
  // (e.g. loops over SlaveStepperList in MasterStepper)
  virtual void clear()
  {
    theOwner->clear();
  }

  virtual void react()
  {
    theOwner->react();
  }

  virtual void turn()
  {
    theOwner->turn();
  }

  virtual void transit()
  {
    theOwner->transit();
  }

  virtual void postern()
  {
    theOwner->postern();
  }

  virtual void initialize()
  {
    Stepper::initialize();
  }

  void setMaster( MasterStepperPtr master ) 
  { 
    theMaster = master; 
  }

  Float getDeltaT() 
  { 
    // Slaves are synchronous to their masters.
    return theMaster->getDeltaT(); 
  }

  static StepperPtr instance() { return new SlaveStepper; }

  virtual const char* const className() const  { return "SlaveStepper"; }

private:

  MasterStepperPtr theMaster;

};


class Euler1Stepper : public MasterStepper
{

public:

  Euler1Stepper();
  virtual ~Euler1Stepper() {}

  static StepperPtr instance() { return new Euler1Stepper; }

  virtual int getNumberOfSteps() { return 1; }
  virtual void clear();
  virtual void react();
//  virtual void check();
  virtual void transit();
  virtual void postern();
  virtual void initialize();

  virtual const char* const className() const  { return "Euler1Stepper"; }
 
protected:

  static IntegratorPtr newEuler1( SubstanceRef substance );

};


class RungeKutta4Stepper : public MasterStepper
{

public:

  RungeKutta4Stepper();
  virtual ~RungeKutta4Stepper() {}

  static StepperPtr instance() { return new RungeKutta4Stepper; }

  virtual int getNumberOfSteps() { return 4; }

  virtual void clear();
  virtual void react();
  virtual void transit();
  virtual void postern();

  virtual void initialize();


  virtual const char* const className() const  { return "RungeKuttaStepper"; }

protected:

  static IntegratorPtr newRungeKutta4( SubstanceRef substance );

};

#endif /* ___STEPPER_H___ */



/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/
