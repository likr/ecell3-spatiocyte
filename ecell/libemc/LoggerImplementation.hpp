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
// written by Masayuki Okayama <smash@e-cell.org> at
// E-CELL Project, Lab. for Bioinformatics, Keio University.
//


#ifndef ___LOGGER_IMPLEMENTATION_H___
#define ___LOGGER_IMPLEMENTATION_H___

#include "libecs/libecs.hpp"
#include "libecs/Logger.hpp" 

namespace libemc
{

  /**
     Pure virtual base class (interface definition) of logger
     implementation.
  */

  class LoggerImplementation
  {

  public:

    LoggerImplementation() {}
    virtual ~LoggerImplementation() {}

    virtual libecs::Logger::DataPointVectorCref
    getData( void ) const = 0;

    virtual libecs::Logger::DataPointVectorCref
    getData( libecs::RealCref start,
	     libecs::RealCref end ) const = 0;

    virtual libecs::Logger::DataPointVectorCref
    getData( libecs::RealCref start,
	     libecs::RealCref end,
	     libecs::RealCref interval ) const = 0;

  };


} // namespace libecs

#endif


