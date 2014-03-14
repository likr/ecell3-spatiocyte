//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 2006-2009 Keio University
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-Cell is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// E-Cell is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with E-Cell -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
// E-Cell Project, Institute for Advanced Biosciences, Keio University.
//

#include <libecs/SpatiocyteThread.hpp>

namespace libecs
{

void Thread::initialize()
{
  runChildren();
  theTotalBoxSize = theStepper.getBoxSize();
  theBoxSize = theTotalBoxSize/theThreadSize;
  theRng.Reseed();
  theMols.resize(theBoxSize);
  theTars.resize(theBoxSize);
  theIDs.resize(theBoxSize);
  theAdjoins.resize(theBoxSize);
  theRands.resize(theBoxSize);
  theAdjBoxes.resize(theBoxSize);
  theAdjAdjBoxes.resize(theBoxSize);
  theAdjMols.resize(theBoxSize * 2 * theTotalBoxSize);
  theAdjTars.resize(theBoxSize * 2 * theTotalBoxSize);
  theAdjAdjMols.resize(theBoxSize);
  theAdjAdjTars.resize(theBoxSize);
  theBorderCounts.resize(theBoxSize * 2 * theTotalBoxSize);
  theBorderMols.resize(theBoxSize * 2 * theTotalBoxSize);
  theBorderTars.resize(theBoxSize * 2 * theTotalBoxSize);
  theRepeatAdjCounts.resize(theBoxSize * theTotalBoxSize);
  theRepeatAdjMols.resize(theBoxSize * theTotalBoxSize);
  theRepeatAdjTars.resize(theBoxSize * theTotalBoxSize);
  for(unsigned i(0); i != theBoxSize; ++i)
    {
      theAdjAdjMols[i].resize(2);
      theAdjAdjTars[i].resize(2);
      for(unsigned j(0); j != 2; ++j)
        {
          theAdjAdjMols[i][j].resize(theTotalBoxSize);
          theAdjAdjTars[i][j].resize(theTotalBoxSize);
        }
    }
  //doWork();
  waitChildren();
}

void Thread::initializeLists()
{
  runChildren();
  for(unsigned i(0); i != theBoxSize; ++i)
    {
      theSpecies[0]->initializeLists((theID*theBoxSize)+i, theRng, theMols[i], theTars[i],
                                     theAdjMols.data() + theTotalBoxSize * 2 * i,
                                     theAdjTars.data() + theTotalBoxSize * 2 * i,
                                     theAdjoins[i], theIDs[i], theAdjBoxes[i],
                                     theAdjAdjBoxes[i], theRands[i]);
    }
  waitChildren();
}

void Thread::doWork()
{
  /*
  for(unsigned i(0); i != theBoxSize; ++i)
    {
      for(unsigned j(0); j != 20000; ++j)
        {
          theBorderMols[0][i].push_back(theRng.IntegerC(11));
          theBorderMols[1][i].push_back(theRng.IntegerC(11));
          theBorderTars[0][i].push_back(theRng.IntegerC(11));
          theBorderTars[1][i].push_back(theRng.IntegerC(11));
        }
    }
    */
}

void Thread::walk()
{
  runChildren();
  unsigned r(0);
  unsigned w(1);
  if(isToggled)
    {
      r = 1;
      w = 0;
      isToggled = false;
    }
  else
    {
      isToggled = true;
    } 
  //if(!theID)
    {
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          theSpecies[0]->walk((theID*theBoxSize)+i, r, w, theRng, theMols[i], theTars[i],
                              theAdjMols.data() + theTotalBoxSize * 2 * i,
                              theAdjTars.data() + theTotalBoxSize * 2 * i,
                              theAdjAdjMols[i], theAdjAdjTars[i],
                              theBorderCounts.data() + theTotalBoxSize * 2 * i,
                              theBorderMols.data() + theTotalBoxSize * 2 * i,
                              theBorderTars.data() + theTotalBoxSize * 2 * i,
                              theRepeatAdjCounts.data() + theTotalBoxSize * i,
                              theRepeatAdjMols.data() + theTotalBoxSize * i,
                              theRepeatAdjTars.data() + theTotalBoxSize * i,
                              theAdjoins[i], theIDs[i], theAdjBoxes[i],
                              theAdjAdjBoxes[i], theRands[i]);
        }
    }
  waitChildren();
}

void Thread::runChildren()
{
  if(!theID)
    {
      if(isRunA)
        {
          flagA = FLAG_RUN;
        }
      else
        {
          flagB = FLAG_RUN;
        }
    }
}

void Thread::waitChildren()
{
  if(!theID)
    {
      while(ACCESS_ONCE(nThreadsRunning) < theThreadSize-1)
        {
          continue;
        }
      nThreadsRunning = 0;
      if(isRunA)
        {
          flagA = FLAG_STOP;
          isRunA = false;
        }
      else
        {
          flagB = FLAG_STOP;
          isRunA = true;
        }
      //__sync_synchronize();
    }
}

void Thread::waitParent()
{
  if(isRunA)
    {
      while(ACCESS_ONCE(flagA) == FLAG_STOP)
        {
          continue;
        }
      isRunA = false;
    }
  else
    {
      while(ACCESS_ONCE(flagB) == FLAG_STOP)
        {
          continue;
        }
      isRunA = true;
    }
}

void Thread::work()
{
  waitParent();
  theStepper.constructLattice(theID);
  __sync_fetch_and_add(&nThreadsRunning, 1);
  waitParent();
  theStepper.concatenateLattice(theID);
  __sync_fetch_and_add(&nThreadsRunning, 1);
  waitParent();
  initialize();
  __sync_fetch_and_add(&nThreadsRunning, 1);
  waitParent();
  initializeLists();
  __sync_fetch_and_add(&nThreadsRunning, 1);
  waitParent();
  /*
  std::vector<std::vector<std::vector<unsigned> > > aBorderMols;
  std::vector<std::vector<std::vector<unsigned> > > aBorderTars;
  aBorderMols.resize(2);
  aBorderTars.resize(2);
  for(unsigned i(0); i != 2; ++i)
    {
      aBorderMols[i].resize(theBoxSize);
      aBorderTars[i].resize(theBoxSize);
    }
  for(unsigned i(0); i != theBoxSize; ++i)
    {
      for(unsigned j(0); j != 20000; ++j)
        {
          aBorderMols[0][i].push_back(theRng.IntegerC(11));
          aBorderMols[1][i].push_back(theRng.IntegerC(11));
          aBorderTars[0][i].push_back(theRng.IntegerC(11));
          aBorderTars[1][i].push_back(theRng.IntegerC(11));
        }
    }
    */
  for(;;)
    {
      //walk(aBorderMols, aBorderTars);
      walk();
      __sync_fetch_and_add(&nThreadsRunning, 1);
      waitParent();
    }
}


void Thread::updateMols(std::vector<unsigned>& aMols, unsigned aBoxID)
{
  aMols.resize(theMols[aBoxID].size());
  for(unsigned i(0); i != theMols[aBoxID].size(); ++i)
    {
      aMols[i] = theMols[aBoxID][i];
    }
  for(unsigned i(0); i != theTotalBoxSize; ++i)
    {
      for(unsigned j(0); j != getAdjMolsSize(i, 0, aBoxID); ++j)
        {
          aMols.push_back(getAdjMols(i, 0, aBoxID)[j]);
        }
      for(unsigned j(0); j != getAdjMolsSize(i, 1, aBoxID); ++j)
        {
          aMols.push_back(getAdjMols(i, 1, aBoxID)[j]);
        }
    }
}

}
