//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-CELL Simulation Environment package
//
//                Copyright (C) 1996-2002 Keio University
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

#include <iostream>

#include "Util.hpp"
#include "VariableReference.hpp"
#include "Stepper.hpp"
#include "FullID.hpp"
#include "Variable.hpp"
#include "Model.hpp"
#include "PropertySlotMaker.hpp"

#include "Process.hpp"


namespace libecs
{

  void Process::makeSlots()
  {
    registerSlot( getPropertySlotMaker()->
		  createPropertySlot( "VariableReferenceList", *this, 
				      Type2Type<Polymorph>(),
				      &Process::setVariableReferenceList,
				      &Process::getVariableReferenceList ) );

    registerSlot( getPropertySlotMaker()->
		  createPropertySlot( "Activity", *this, 
				      Type2Type<Real>(),
				      &Process::setActivity,
				      &Process::getActivity ) );

    registerSlot( getPropertySlotMaker()->
		  createPropertySlot( "Priority", *this, 
				      Type2Type<Int>(),
				      &Process::setPriority,
				      &Process::getPriority ) );
  }


  void Process::setVariableReferenceList( PolymorphCref aValue )
  {
    const PolymorphVector aVector( aValue.asPolymorphVector() );
    for( PolymorphVectorConstIterator i( aVector.begin() );
	 i != aVector.end(); ++i )
      {
	const PolymorphVector anInnerVector( (*i).asPolymorphVector() );

	setVariableReference( anInnerVector );
      }

  }

  const Polymorph Process::getVariableReferenceList() const
  {
    PolymorphVector aVector;
    aVector.reserve( theVariableReferenceVector.size() );
  
    for( VariableReferenceVectorConstIterator i( theVariableReferenceVector.begin() );
	 i != theVariableReferenceVector.end() ; ++i )
      {
	PolymorphVector anInnerVector;
	VariableReferenceCref aVariableReference( *i );

	// Tagname
	anInnerVector.push_back( aVariableReference.getName() );
	// FullID
	anInnerVector.push_back( aVariableReference.getVariable()->
				 getFullID().getString() );
	// Coefficient
	anInnerVector.push_back( aVariableReference.getCoefficient() );

	aVector.push_back( anInnerVector );
      }

    return aVector;
  }


  Process::Process() 
    :
    theFirstZeroVariableReference( theVariableReferenceVector.end() ),
    theFirstPositiveVariableReference( theVariableReferenceVector.end() ),
    theActivity( 0.0 ),
    thePriority( 0 )
  {
    makeSlots();
  }

  Process::~Process()
  {
    ; // do nothing
  }


  VariableReference Process::getVariableReference( StringCref aName )
  {
    return *( findVariableReference( aName ) );
  }

  void Process::removeVariableReference( StringCref aName )
  {
    theVariableReferenceVector.erase( findVariableReference( aName ) );
  }

  void Process::setVariableReference( PolymorphVectorCref aValue )
  {

    UnsignedInt aVectorSize( aValue.size() );
    
    // Require ( tagname, fullid, coefficient ) 3-tuple
    if( aVectorSize < 2 )
      {
	THROW_EXCEPTION( ValueError, "Process [" + getFullID().getString()
			 + "]: ill-formed VariableReference given." );
      }

    const String aVariableReferenceName(  aValue[0].asString() );
    const String aFullIDString( aValue[1].asString() );
    if( ! aFullIDString.empty() )
      {
	const FullID aFullID( aValue[1].asString() );
	Int          aCoefficient( 0 );
	
	if( aVectorSize >= 3 )
	  {
	    aCoefficient = aValue[2].asInt();
	  }
	
	registerVariableReference( aVariableReferenceName, aFullID, 
				   aCoefficient );
      }
    else // if the FullID is empty, remove the VariableReference
      {
	removeVariableReference( aVariableReferenceName );
      }
  }



  void Process::registerVariableReference( StringCref aName, 
					   FullIDCref aFullID, 
					   const Int aCoefficient )
  {
    SystemPtr aSystem( getModel()->getSystem( aFullID.getSystemPath() ) );
    VariablePtr aVariable( aSystem->getVariable( aFullID.getID() ) );

    registerVariableReference( aName, aVariable, aCoefficient );
  }


  void Process::registerVariableReference( StringCref aName, 
					   VariablePtr aVariable, 
					   const Int aCoefficient )
  {
    VariableReference aVariableReference( aName, aVariable, aCoefficient );
    theVariableReferenceVector.push_back( aVariableReference );

    // sort by coefficient
    std::sort( theVariableReferenceVector.begin(), 
	       theVariableReferenceVector.end(), 
	       VariableReference::CoefficientCompare() );

    // find the first VariableReference whose coefficient is 0,
    // and the first VariableReference whose coefficient is positive.
    std::pair<VariableReferenceVectorConstIterator,
      VariableReferenceVectorConstIterator> 
      aZeroRange( std::equal_range( theVariableReferenceVector.begin(), 
				    theVariableReferenceVector.end(), 
				    0, 
				    VariableReference::CoefficientCompare()
				    ) );

    theFirstZeroVariableReference     = aZeroRange.first;
    theFirstPositiveVariableReference = aZeroRange.second;
  }


  VariableReferenceVectorIterator 
  Process::findVariableReference( StringCref aName )
  {
    // well this is a linear search.. but this won't be used in simulation.
    for( VariableReferenceVectorIterator 
	   i( theVariableReferenceVector.begin() );
	 i != theVariableReferenceVector.end(); ++i )
      {
	if( (*i).getName() == aName )
	  {
	    return i;
	  }
      }

    THROW_EXCEPTION( NotFound,
		     "[" + getFullID().getString() + 
		     "]: VariableReference [" + aName + 
		     "] not found in this Process." );
  }


  void Process::initialize()
  {
    ; // do nothing
  }


} // namespace libecs


/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/
