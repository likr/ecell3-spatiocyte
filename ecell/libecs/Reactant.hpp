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

#ifndef __REACTANT_H___
#define __REACTANT_H___
#include "Substance.hpp"
#include "libecs.hpp"

class Reactant
{

public:

  Reactant( SubstanceRef s, const int c ) 
    : 
    theSubstance( s ), theCoefficient( c ) 
    {
      ; // do nothing
    }
  virtual ~Reactant() {}

  SubstanceRef getSubstance() const { return theSubstance; }
  int getCoefficient() const { return theCoefficient; }
  Float getConcentration() const 
  { return theSubstance.getConcentration(); }
  Float getQuantity() const { return theSubstance.getQuantity(); }
  Float getActivity() const { return theSubstance.getActivity(); }
  Float getVelocity() const { return theSubstance.getVelocity(); }
  void addVelocity( Float v ) const { theSubstance.addVelocity( v ); }
  void setQuantity( Float q ) const { theSubstance.setQuantity( q ); }

private:

  SubstanceRef theSubstance;
  int theCoefficient;

};

#endif /* __REACTANT_H___ */
