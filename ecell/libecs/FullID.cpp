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

#include <string>

#include "Util.hpp"
#include "Exceptions.hpp"

#include "FullID.hpp"

namespace libecs
{

  ///////////////////////  SystemPath

  void SystemPath::parse( StringCref systempathstring )
  {
    if( systempathstring == "" )
      {
	return;
      }

    String aString( systempathstring );
    eraseWhiteSpaces( aString );

    // ignore leading white spaces
    String::size_type aFieldStart( 0 );
    String::size_type aFieldEnd( aString.
				 find_first_of( DELIMITER,
						aFieldStart ) );
    
     // absolute path ( start with '/' )
    if( aFieldEnd == 0 )
       {
	 push_back( String( 1, DELIMITER ) );
	 ++aFieldStart;
	 aFieldEnd = aString.find_first_of( DELIMITER, aFieldStart );

	 if( aFieldEnd == String::npos )
	   {
	     return;
	   }
       }

    push_back( aString.substr( aFieldStart, 
			       aFieldEnd - aFieldStart ) );

    while( aFieldEnd != String::npos  )
      {
	aFieldStart = aFieldEnd + 1;
	aFieldEnd = aString.find_first_of( DELIMITER, aFieldStart );

	push_back( aString.substr( aFieldStart, 
				   aFieldEnd - aFieldStart ) );
      }

  }

  const String SystemPath::getString() const
  {
    StringListConstIterator i = begin();
    String aString;

    if( isAbsolute() )
      {
	if( size() == 1 )
	  {
	    return "/";
	  }
      }
    else
      {
	// isAbsolute() == false implies that it can be empty
	if( i == end() )
	  {
	    return "";
	  }

	aString = *i;
      }


    ++i;

    while( i != end() )
      {
	aString += '/';
	aString += *i;
	++i;
      }

    return aString;
  }



  ///////////////// FullID

  void FullID::parse( StringCref fullidstring )
  {
    // empty FullID string is invalid
    if( fullidstring == "" )
      {
	throw BadID( __PRETTY_FUNCTION__,
		     "Empty FullID string." );
      }

    String aString( fullidstring );
    eraseWhiteSpaces( aString );

    // ignore leading white spaces
    String::size_type aFieldStart( 0 );
    String::size_type aFieldEnd( aString.find_first_of( DELIMITER,
							     aFieldStart ) );
    if( aFieldEnd == String::npos )
      {
	throw BadID( __PRETTY_FUNCTION__,
		     "No ':' in the FullID string [" + aString + "]." );
      }

    String aTypeString( aString.substr( aFieldStart, 
					     aFieldEnd - aFieldStart ) );
    thePrimitiveType = PrimitiveTypeOf( aTypeString );

    aFieldStart = aFieldEnd + 1;
    aFieldEnd = aString.find_first_of( DELIMITER, aFieldStart );
    if( aFieldEnd == String::npos )
      {
	throw BadID( __PRETTY_FUNCTION__,
		     "Only one ':' in the FullID string [" 
		     + aString + "]." );
      }

    theSystemPath = 
      SystemPath( aString.substr( aFieldStart, 
				       aFieldEnd - aFieldStart ) );

    aFieldStart = aFieldEnd + 1;

    // drop trailing string after extra ':'(if this is  FullPropertyName),
    // or go to the end
    aFieldEnd = aString.find_first_of( DELIMITER, aFieldStart );

    theID = aString.substr( aFieldStart, aFieldEnd - aFieldStart );
  }    

  const String FullID::getString() const
  {
    return PrimitiveTypeStringOf( thePrimitiveType ) + FullID::DELIMITER 
      + theSystemPath.getString() + FullID::DELIMITER + theID;
  }

  bool FullID::isValid() const
  {
    bool aFlag( theSystemPath.isValid() );
    aFlag &= ! theID.empty();

    return aFlag;
  }


  ///////////////// FullPropertyName


  FullPropertyName::FullPropertyName( StringCref fullpropertynamestring )
    :
    theFullID( fullpropertynamestring )
  {

    String::size_type aPosition( 0 );

    for( Int i( 0 ) ; i < 3 ; ++i )
      {
	aPosition = fullpropertynamestring.
	  find_first_of( FullID::DELIMITER, aPosition );
	if( aPosition == String::npos ) 
	  {
	    throw BadID( __PRETTY_FUNCTION__, 
			 "Not enough fields in FullPropertyName string [" +
			 fullpropertynamestring + "]." );
	  }
	++aPosition;
      }

    thePropertyName = fullpropertynamestring.substr( aPosition, String::npos );
    eraseWhiteSpaces( thePropertyName );
  }

  const String FullPropertyName::getString() const
  {
    return theFullID.getString() + FullID::DELIMITER + thePropertyName;
  }

  bool FullPropertyName::isValid() const
  {
    return theFullID.isValid() & ! thePropertyName.empty();
  }

} // namespace libecs

#ifdef TEST_FQPI

using namespace libecs;

main()
{
  SystemPath aSystemPath( "   \t  /A/BB/CCC//DDDD/EEEEEE    \t \n  " );
  cout << aSystemPath.getString() << endl;

  SystemPath aSystemPath2( aSystemPath );
  cout << aSystemPath2.getString() << endl;
  
  aSystemPath2.pop_front();
  aSystemPath2.pop_back();
  cout << aSystemPath2.getString() << endl;

  SystemPath aSystemPath3( "/" );
  cout << aSystemPath3.getString() << endl;
  cout << aSystemPath3.size() << endl;
  aSystemPath3.pop_front();
  cout << aSystemPath3.getString() << endl;
  cout << aSystemPath3.size() << endl;
  cout << aSystemPath3.empty() << endl;

  while( aSystemPath.size() != 0 )
    {
      cout << aSystemPath.size() << " : " << aSystemPath.isAbsolute() << " : " 
	   << aSystemPath.getString() << endl;
      aSystemPath.pop_front();
    }

  //  SystemPath aSystemPath2( "/A/../B" );
  //  cout << aSystemPath2.getString() << endl;

  cout << "\n::::::::::" << endl;

  try
    {
      FullID aFullID( "       \t  \n  Substance:/A/B:S   \t   \n" );
      cout << aFullID.getString() << endl;
      cout << aFullID.getPrimitiveType() << endl;
      cout << aFullID.getSystemPath().getString() << endl;
      cout << aFullID.getID() << endl;
      cout << aFullID.isValid() << endl;

      FullID aFullID2( aFullID );
      cout << aFullID2.getString() << endl;

      FullID aFullID3( "Reactor:/:R" );
      cout << aFullID3.getString() << endl;
      aFullID3 = aFullID2;
      cout << aFullID3.getString() << endl;

      cout << "\n::::::::::" << endl;

      FullPropertyName 
	aFullPropertyName( "       \t  \n  Substance:/A/B:S:PNAME   \t   \n" );
      cout << aFullPropertyName.getString() << endl;
      cout << aFullPropertyName.getPrimitiveType() << endl;
      cout << aFullPropertyName.getSystemPath().getString() << endl;
      cout << aFullPropertyName.getID() << endl;
      cout << aFullPropertyName.getPropertyName() << endl;
      cout << aFullPropertyName.isValid() << endl;

      FullPropertyName aFullPropertyName2( aFullPropertyName );
      cout << aFullPropertyName2.getString() << endl;

      FullPropertyName aFullPropertyName3( "Reactor:/:R:P" );
      cout << aFullPropertyName3.getString() << endl;
      aFullPropertyName3 = aFullPropertyName2;
      cout << aFullPropertyName3.getString() << endl;

    }
  catch ( ExceptionCref e )
    {
      cerr << e.message() << endl;
    }

}


#endif


/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/
