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

#include "Integrators.hpp"

////////////////////////////// Integrator

Integrator::Integrator( SubstanceRef substance ) 
  :
  theSubstance( substance ),
  theStepCounter( 0 )
{
  theSubstance.setIntegrator( this );
}

void Integrator::clear()
{
  theStepCounter=0;
  theNumberOfMoleculesCache = theSubstance.getQuantity();
}

////////////////////////////// Eular1Integrator

Eular1Integrator::Eular1Integrator( SubstanceRef substance ) 
  : 
  Integrator( substance )
{
  ; // do nothing
}

void Eular1Integrator::turn()
{
  ++theStepCounter;
}


////////////////////////////// RungeKutta4Integrator

const Float RungeKutta4Integrator::theOne6th = (Float)1.0 / (Float)6.0;

RungeKutta4Integrator::TurnFunc RungeKutta4Integrator::theTurnFuncs[4] =
{
  &RungeKutta4Integrator::turn0,
  &RungeKutta4Integrator::turn1,
  &RungeKutta4Integrator::turn2,
  &RungeKutta4Integrator::turn3
};


RungeKutta4Integrator::RungeKutta4Integrator( SubstanceRef substance ) 
  : 
  Integrator(substance)
{
  ; // do nothing
}

void RungeKutta4Integrator::clear()
{
  Integrator::clear();
  theTurnFuncPtr = &theTurnFuncs[0];
}

void RungeKutta4Integrator::turn()
{
  theK[ theStepCounter ] = theSubstance.getVelocity();
  ( this->*( *theTurnFuncPtr ) )();
  ++theTurnFuncPtr;
  ++theStepCounter;
  setVelocity( 0 );
}

void RungeKutta4Integrator::turn0()
{
  setQuantity( ( theK[0] * .5 ) + theNumberOfMoleculesCache );
}

void RungeKutta4Integrator::turn1()
{
  setQuantity( ( theK[1] * .5 ) + theNumberOfMoleculesCache );
}

void RungeKutta4Integrator::turn2()
{
  setQuantity( theK[2] + theNumberOfMoleculesCache );
}

void RungeKutta4Integrator::turn3()
{
  setQuantity( theNumberOfMoleculesCache );
}

void RungeKutta4Integrator::transit()
{
  //// x(n+1) = x(n) + 1/6 * (k1 + k4 + 2 * (k2 + k3)) + O(h^5)

  Float* k( &theK[0] );

  // FIXME: prefetching here makes this faster on alpha?

  //                      k1  + (2 * k2)  + (2 * k3)  + k4
  Float aResult = *k++;
  aResult += *k;
  aResult += *k++;
  aResult += *k;
  aResult += *k++;
  aResult += *k;

  aResult *= theOne6th;

  setVelocity( aResult );
}


/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/
