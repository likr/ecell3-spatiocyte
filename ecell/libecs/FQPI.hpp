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

#ifndef ___FQPI_H___
#define ___FQPI_H___
#include <string>

#include "libecs.hpp"
#include "Exceptions.hpp"
#include "PrimitiveType.hpp"

/** 
  SystemPath 
  */
class SystemPath {

public:

  static const char DELIMITER = '/';

  // exceptions.

  class SystemPathException : public Exception
    { 
    public: 
      SystemPathException( StringCref method, StringCref what ) 
	: Exception( method, what ) {} 
      const String what() const { return ""; }
    };
  class BadSystemPath : public SystemPathException
    { 
    public: 
      BadSystemPath( StringCref method, StringCref what ) 
	: SystemPathException( method, what ) {} 
      const String what() const { return "Bad SystemPath."; }
    };

public:

  SystemPath( StringCref systempath = "" );
  virtual ~SystemPath() {}

  StringCref getSystemPath() const { return theSystemPath; }
  virtual const String getString() const { return getSystemPath(); }

  virtual operator String() const { return getString(); }

  /**
    Extract the first system name. Standardize given string.
    @return name of the first system
    */
  const String first() const;

  /**
    Extract the last system name. Standardize given string.

    @return name of the last system in given systempath.
    */
  const String last() const;
  /**
    Remove the first system name. Standardize given string.
    @return
    */
  SystemPath next() const;

protected:

  /**
    Standardize a SystemPath. 
    Reduce '..'s and remove extra white spaces.

    @return reference to the systempath
    */
  void standardize();

  SystemPath() {}

private:

  const String theSystemPath;

};

/**
  FQID(Fully Qualified entity ID)

  The FQID is a identifier (ID) of Entity objects of certain Primitive
  type.  Given a Primitive type, one can identify unique Entity in a
  cell model with a SystemPath and an id.  

  @see SystemPath, Primitive 
*/
class FQID : public SystemPath
{
public: // exceptions

  class FQIDException : public Exception
    { 
    public: 
      FQIDException( StringCref method, StringCref what )
	: Exception( method, what ) {} 
      const String what() const { return ""; }
    };

  class BadFQID : public FQIDException
    { 
    public: 
      BadFQID( StringCref method, StringCref what ) 
	: FQIDException( method,what ) {} 
      const String what() const { return "Bad FQID"; }
    };

public:

  FQID( StringCref systemname, StringCref id );
  FQID( StringCref fqen );
  virtual ~FQID() {}

  const String getFqid() const;
  virtual const String getString() const { return getFqid(); }
  StringCref getId() const { return theId; }
  virtual operator String() const { return getString(); }

  static const String IdOf( StringCref fqen );
  static const String SystemPathOf( StringCref fqen );

private:

  const String theId;

};

class FQPIException : public Exception
{ 
public: 
  FQPIException( StringCref method, StringCref message ) 
    : Exception( method, message ) {} 
  const String what() const { return ""; }
};
class BadFQPI : public FQPIException
{ 
public:
  BadFQPI( StringCref method, StringCref message ) 
    : FQPIException( method, message ) {} 
  const String what() const { return "Bad FQPI."; }
};



/**
  FQPI (Fully Qualified Primitive Id).

  One can identify an unique Entiy in a cell model with a FQPI.
  The FQPI consists of FQID and PrimitiveType.

  @see FQID, PrimitiveType
*/
class FQPI : public FQID
{

public:

  static const String  fqidOf( StringCref fqpi );


  FQPI( const PrimitiveType type, FQIDCref fqid );
  FQPI( StringCref fqpi );
  virtual ~FQPI() {}
  
  const String getFqpi() const;
  const PrimitiveType getPrimitiveType() const { return thePrimitiveType; }

  virtual const String getString() const { return getFqpi(); }
  virtual operator const String() const { return getString(); }

private:

  PrimitiveType thePrimitiveType;

};

#endif /*  ___FQPI_H___ */

/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/
