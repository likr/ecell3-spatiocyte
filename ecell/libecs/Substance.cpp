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

#include "Substance.hpp"
#include "Integrators.hpp"
#include "System.hpp"
#include "RootSystem.hpp"
#include "Accumulators.hpp"
#include "AccumulatorMaker.hpp"
#include "PrimitiveType.hpp"
#include "Util.hpp"

const String Substance::SYSTEM_DEFAULT_ACCUMULATOR_NAME = "ReserveAccumulator";

String Substance::USER_DEFAULT_ACCUMULATOR_NAME 
= Substance::SYSTEM_DEFAULT_ACCUMULATOR_NAME;

void Substance::makeSlots()
{
  MessageSlot( "id",Substance,*this,NULL,
	       &Substance::getId );
  MessageSlot( "quantity",Substance,*this,&Substance::setQuantity,
	      &Substance::getQuantity );
  MessageSlot( "accumulator",Substance,*this,&Substance::setAccumulator,
	       &Substance::getAccumulator );
}

void Substance::setQuantity( MessageCref message )
{
  Float aQuantity = stringTo<Float>( message.getBody() );

  if( theAccumulator )
    {
      loadQuantity( aQuantity );
    }
  else
    {
      mySetQuantity( aQuantity );
    }
}


void Substance::setAccumulator( MessageCref message )
{
  setAccumulator( message.getBody() );
}

const Message Substance::getQuantity( StringCref keyword )
{
  static String aKeyword( "quantity" );
  return Message( aKeyword, getQuantity() );
}

const Message Substance::getAccumulator( StringCref keyword )
{
  static String aKeyword( "accumulator" );
  if( theAccumulator )
    {
      return Message( aKeyword, theAccumulator->className() );
    }
  else
    {
      return Message( aKeyword, "" );
    }
}

Substance::Substance()
  : 
  theAccumulator( NULL ),
  theIntegrator( NULL ),
  theQuantity( 0 ),  
  theFraction( 0 ),
  theVelocity( 0 ),
  theBias( 0 ),
  theFixed( false ) ,
  theConcentration( -1 )
{
  makeSlots();
} 

Substance::~Substance()
{
  delete theIntegrator;
  delete theAccumulator;
}

void Substance::setAccumulator( StringCref classname )
{
  try {
    AccumulatorPtr aAccumulatorPtr;
    aAccumulatorPtr = getSuperSystem()->getRootSystem()->
      getAccumulatorMaker().make( classname );
    setAccumulator(aAccumulatorPtr);
    if( classname != userDefaultAccumulatorName() )
      {
	cerr << __PRETTY_FUNCTION__ << endl;
	//FIXME:    *theMessageWindow << "[" << fqpi() 
	//FIXME: << "]: accumulator is changed to: " << classname << ".\n";
      }
  }
  catch( Exception& e )
    {
      //FIXME:     *theMessageWindow << "[" << fqpi() << "]:\n" << e.message();
      // warn if theAccumulator is already set
      if( theAccumulator != NULL )   
       {
	 //FIXME: *theMessageWindow << "[" << fqpi() << 
	 //FIXME: "]: falling back to :" << theAccumulator->className() 
	 //FIXME: << ".\n";
       }
  }
}

void Substance::setAccumulator( AccumulatorPtr accumulator )
{
  if( theAccumulator )
    delete theAccumulator;
  theAccumulator = accumulator;
  theAccumulator->setOwner( this );
  theAccumulator->update();
}

const String Substance::getFqpi() const
{
  return PrimitiveTypeStringOf( *this ) + ":" + getFqid();
}


void Substance::initialize()
{
  if( theAccumulator == NULL )
    setAccumulator( USER_DEFAULT_ACCUMULATOR_NAME );

  // if the user default is invalid fall back to the system default.
  if( !theAccumulator )  
    {               
      //FIXME:      *theMessageWindow << "Substance: " 
      //FIXME:	<< "falling back to the system default accumulator: " 
      //FIXME:	  << SYSTEM_DEFAULT_ACCUMULATOR_NAME  << ".\n";
      setUserDefaultAccumulatorName(SYSTEM_DEFAULT_ACCUMULATOR_NAME);
      setAccumulator(USER_DEFAULT_ACCUMULATOR_NAME);
    }
}

Float Substance::saveQuantity()
{
  return theAccumulator->save();
}

void Substance::loadQuantity( Float quantity )
{
  mySetQuantity( quantity );
  theAccumulator->update();
}

Float Substance::getActivity()
{
  return getVelocity();
}

void Substance::calculateConcentration()
{
  theConcentration = theQuantity / ( getSuperSystem()->getVolume() * N_A ); 
}


bool Substance::haveConcentration() const
{
  bool aBool(true);

  if( getSuperSystem()->getVolumeIndex() == NULL ) 
    {
      aBool = false;
    }

  return aBool;
}

void Substance::transit()
{ 
  theConcentration = -1;

  if( theFixed ) 
    return;

  theIntegrator->transit();

  theAccumulator->doit();
  
   if( theQuantity < 0 ) 
     {
       theQuantity = 0;
//       throw LTZ();
     }
}

void Substance::clear()
{ 
  theQuantity += theBias; 
  theBias = 0;
  theVelocity = 0; 
  theIntegrator->clear();
}

void Substance::turn()
{
  theIntegrator->turn();
}


/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/
