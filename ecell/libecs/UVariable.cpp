//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-CELL Simulation Environment package
//
//                Copyright (C) 1996-2001 Keio University
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

#include "Util.hpp"

#include "UniversalVariable.hpp"


UniversalVariableStringData::
UniversalVariableStringData( const Float f )
  :
  theString( toString<Float>( f ) )
{
  ; // do nothing
}

UniversalVariableStringData::UniversalVariableStringData( const Int i )
  :
  theString( toString<Int>( i ) )
{
  ; // do nothing
}

const Float UniversalVariableStringData::asFloat() const
{
  return stringTo<Float>( theString );
}

const Int UniversalVariableStringData::asInt() const
{
  return stringTo<Int>( theString );
}


UniversalVariableFloatData::
UniversalVariableFloatData( StringCref str )
  :
  theFloat( stringTo<Float>( str ) )
{
  ; // do nothing
}

const String UniversalVariableFloatData::asString() const
{
  return toString<Float>( theFloat );
}

UniversalVariableIntData::
UniversalVariableIntData( StringCref str )
  :
  theInt( stringTo<Int>( str ) )
{
  ; // do nothing
}

UniversalVariableIntData::
UniversalVariableIntData( const Float f )
  :
  // FIXME: range check?
  theInt( static_cast<Int>( f ) )
{
  ; // do nothing
}


const String UniversalVariableIntData::asString() const
{
  return toString<Int>( theInt );
}


