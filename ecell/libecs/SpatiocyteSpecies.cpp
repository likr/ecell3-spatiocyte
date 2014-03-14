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

#include <SpatiocyteSpecies.hpp>
#if defined(__MIC__)
#include <immintrin.h>
#endif

namespace libecs
{

//all writing -> local
//all reading -> shared

void Species::updateBoxMols(const unsigned currBox, const unsigned r,
                            std::vector<unsigned>& aMols,
                            std::vector<unsigned>& aTars,
                            const std::vector<unsigned>& anAdjBoxes)
{
  const unsigned aBoxSizePerThread = theBoxSize / theThreadSize;
  for(unsigned i(0); i != anAdjBoxes.size(); ++i)
    {
      const unsigned adjBox(anAdjBoxes[i] / aBoxSizePerThread);
      const unsigned aBoxID(anAdjBoxes[i] % aBoxSizePerThread);
      Thread* aThread(theThreads[adjBox]);
      //reading border mols, so shared:
      if (aThread->getBorderCount(currBox, r, aBoxID))
        {
          std::vector<unsigned>& borderMols(aThread->getBorderMols(currBox, r, aBoxID));
          std::vector<unsigned>& borderTars(aThread->getBorderTars(currBox, r, aBoxID));
          for(unsigned j(0); j < borderMols.size(); ++j)
            {
              aMols.push_back(borderMols[j]);
              aTars.push_back(borderTars[j]);
            }
          aThread->setBorderCount(currBox, r, aBoxID, 0);
          borderMols.resize(0);
          borderTars.resize(0);
        }
    }
}

void Species::walkMols(std::vector<unsigned>& aMols,
                       const std::vector<unsigned>& aTars,
                       std::vector<unsigned short>& anIDs)
{
#if defined(__MIC__)
  const unsigned BLOCKSIZE(64);
  const unsigned aSize(aMols.size());
  for(unsigned iStart(0); iStart < aSize; iStart += BLOCKSIZE)
    {
      const unsigned iStop = std::min(iStart + BLOCKSIZE, aSize);
      for(unsigned i(iStart); i < iStop; ++i)
        {
          const unsigned aTar(aTars[i]);
          const unsigned aTarMol(aTar%theBoxMaxSize);
          _mm_prefetch((const char*)(anIDs.data() + aTarMol), _MM_HINT_T0);
          //_mm_prefetch((const char*)(anIDs.data() + aMols[i]), _MM_HINT_T0);
        }
      for(unsigned i(iStart); i < iStop; ++i)
        {
          const unsigned aTar(aTars[i]);
          const unsigned aTarMol(aTar%theBoxMaxSize);
          if(anIDs[aTarMol] == theVacantID)
            {
              anIDs[aTarMol] = theID;
              anIDs[aMols[i]] = theVacantID;
              aMols[i] = aTarMol;
            }
        }
    }
#else
  for(unsigned i(0); i < aMols.size(); ++i)
    {
      const unsigned aTar(aTars[i]);
      const unsigned aTarMol(aTar%theBoxMaxSize);
      if(anIDs[aTarMol] == theVacantID)
        {
          anIDs[aTarMol] = theID;
          anIDs[aMols[i]] = theVacantID;
          aMols[i] = aTarMol;
        }
    }
#endif
}

void Species::updateAdjMols(const unsigned currBox, const unsigned r,
                            std::vector<std::vector<unsigned> >& aRepeatAdjMols,
                            std::vector<std::vector<unsigned> >& aRepeatAdjTars,
                            const std::vector<unsigned>& anAdjBoxes)
{
  for(unsigned i(0); i != anAdjBoxes.size(); ++i)
    {
      const unsigned adjBox(anAdjBoxes[i]/(theBoxSize/theThreadSize));
      const unsigned aBoxID(anAdjBoxes[i]%(theBoxSize/theThreadSize));
      std::vector<unsigned>& repeatAdjMols(aRepeatAdjMols[anAdjBoxes[i]]);
      std::vector<unsigned>& repeatAdjTars(aRepeatAdjTars[anAdjBoxes[i]]);
      std::vector<unsigned>& adjMols(theThreads[adjBox
                                     ]->getAdjMols(currBox, r, aBoxID));
      std::vector<unsigned>& adjTars(theThreads[adjBox
                                     ]->getAdjTars(currBox, r, aBoxID));
      for(unsigned j(0); j != repeatAdjMols.size(); ++j)
        {
          adjMols.push_back(repeatAdjMols[j]);
          adjTars.push_back(repeatAdjTars[j]);
        }
      repeatAdjMols.resize(0);
      repeatAdjTars.resize(0);
    }
}

void Species::updateAdjAdjMols(const unsigned currBox, const unsigned r,
                               const std::vector<unsigned>& anAdjAdjBoxes)
{
  //for(unsigned i(0); i != theBoxSize; ++i)
  for(unsigned i(0); i != anAdjAdjBoxes.size(); ++i)
    {
      const unsigned adjAdjBox(anAdjAdjBoxes[i]/(theBoxSize/theThreadSize));
      const unsigned aBoxID(anAdjAdjBoxes[i]%(theBoxSize/theThreadSize));
      //const unsigned adjAdjBox(i);
      std::vector<unsigned>& adjAdjMols(theThreads[adjAdjBox
                                        ]->getAdjAdjMols(currBox, r, aBoxID));
      std::vector<unsigned>& adjAdjTars(theThreads[adjAdjBox
                                        ]->getAdjAdjTars(currBox, r, aBoxID));
      for(unsigned j(0); j != adjAdjMols.size(); ++j)
        {
          const unsigned aMol(adjAdjMols[j]);
          const unsigned aBox(aMol/theBoxMaxSize);
          theThreads[aBox/(theBoxSize/theThreadSize)]->pushAdj(currBox, r,
       aMol-theBoxMaxSize*aBox, adjAdjTars[j], aBox%(theBoxSize/theThreadSize));
        }
      adjAdjMols.resize(0);
      adjAdjTars.resize(0);
    }
}

void Species::walkAdjMols(const unsigned currBox, const unsigned r,
                          std::vector<unsigned>& aMols,
                          std::vector<unsigned short>& anIDs,
                          const std::vector<unsigned>& anAdjBoxes)
{
#if defined(__MIC__)
  const unsigned BLOCKSIZE(512);
  const unsigned aBoxSizePerThread = theBoxSize / theThreadSize;
  for(unsigned iStart(0); iStart < anAdjBoxes.size(); iStart += BLOCKSIZE)
    {
      const unsigned iStop = std::min(iStart + BLOCKSIZE, (unsigned)anAdjBoxes.size());
      for(unsigned i(iStart); i < iStop; ++i)
        {
          const unsigned anAdjBox = anAdjBoxes[i];
          const unsigned aBox(anAdjBox / aBoxSizePerThread);
          const unsigned aBoxID(anAdjBox % aBoxSizePerThread);
          _mm_prefetch((const char*)theThreads[aBox]->getAdjMolsAddress(currBox, r, aBoxID), _MM_HINT_T2);
          _mm_prefetch((const char*)theThreads[aBox]->getAdjTarsAddress(currBox, r, aBoxID), _MM_HINT_T2);
        }
      for(unsigned i(iStart); i < iStop; ++i)
        {
          const unsigned anAdjBox = anAdjBoxes[i];
          const unsigned aBox(anAdjBox / aBoxSizePerThread);
          const unsigned aBoxID(anAdjBox % aBoxSizePerThread);
          std::vector<unsigned>& adjMols(theThreads[aBox]->getAdjMols(currBox, r, aBoxID));
          std::vector<unsigned>& adjTars(theThreads[aBox]->getAdjTars(currBox, r, aBoxID));
          unsigned n(adjMols.size());
          for(unsigned j(0); j < n; ++j)
            {
              const unsigned aTar(adjTars[j]);
              const unsigned aTarMol(aTar%theBoxMaxSize);
              if(anIDs[aTarMol] == theVacantID)
                {
                  anIDs[aTarMol] = theID;
                  theThreads[aBox]->setMolID(adjMols[j], theVacantID, aBoxID);
                  aMols.push_back(aTarMol);
                  adjMols[j] = adjMols.back();
                  adjMols.pop_back();
                  adjTars[j] = adjTars.back();
                  adjTars.pop_back();
                  --j;
                }
            }
          adjTars.resize(0);
        }
    }
#else
  for(unsigned i(0); i != anAdjBoxes.size(); ++i)
    {
      const unsigned anAdjBox = anAdjBoxes[i];
      const unsigned aBox(anAdjBox/(theBoxSize/theThreadSize));
      const unsigned aBoxID(anAdjBox%(theBoxSize/theThreadSize));
      std::vector<unsigned>& adjMols(theThreads[aBox]->getAdjMols(currBox, r,
                                                                  aBoxID));
      std::vector<unsigned>& adjTars(theThreads[aBox]->getAdjTars(currBox, r,
                                                                  aBoxID));
      for(unsigned j(0), n(adjMols.size()); j < n; ++j)
        {
          const unsigned aTar(adjTars[j]);
          const unsigned aTarMol(aTar%theBoxMaxSize);
          if(anIDs[aTarMol] == theVacantID)
            {
              anIDs[aTarMol] = theID;
              theThreads[aBox]->setMolID(adjMols[j], theVacantID, aBoxID);
              //theIDs[aBox][adjMols[j]] = theVacantID;
              aMols.push_back(aTarMol);
              adjMols[j] = adjMols.back();
              adjMols.pop_back();
              adjTars[j] = adjTars.back();
              adjTars.pop_back();
              --j;
            }
        }
      adjTars.resize(0);
    }
#endif
}

void Species::setAdjTars(const unsigned currBox, const unsigned r,
                unsigned* aBorderCounts,
                std::vector<unsigned>* aBorderMols,
                std::vector<unsigned>* aBorderTars,
                std::vector<std::vector<unsigned> >& anAdjAdjMols,
                std::vector<std::vector<unsigned> >& anAdjAdjTars,
                std::vector<std::vector<unsigned> >& aRepeatAdjMols,
                std::vector<std::vector<unsigned> >& aRepeatAdjTars,
                const std::vector<unsigned>& anAdjBoxes,
                unsigned startRand,
                std::vector<unsigned>& aRands)

{
  for(unsigned i(0); i != anAdjBoxes.size(); ++i)
    {
      const unsigned adjBox(anAdjBoxes[i]);
      const unsigned aBoxID(anAdjBoxes[i]%(theBoxSize/theThreadSize));
      //reading adjMols, so get it from the thread:
      std::vector<unsigned>& adjMols(theThreads[anAdjBoxes[i]/
               (theBoxSize/theThreadSize)]->getAdjMols(currBox, r, aBoxID));
      const std::vector<unsigned>& anAdjoins(theAdjoins[adjBox]);
      for(unsigned j(0); j < adjMols.size(); ++j)
        {
          const unsigned aMol(adjMols[j]);
          const unsigned aTar(anAdjoins[aMol*theAdjoinSize+aRands[startRand++]]);
          const unsigned aBox(aTar/theBoxMaxSize);
          if(aBox == currBox)
            {
              aRepeatAdjMols[adjBox].push_back(aMol);
              aRepeatAdjTars[adjBox].push_back(aTar);
            }
          else if(aBox == adjBox)
            {
              aBorderCounts[adjBox] += 1;
              aBorderMols[adjBox].push_back(aMol);
              aBorderTars[adjBox].push_back(aTar);
            }
          else
            {
              anAdjAdjMols[aBox].push_back(theBoxMaxSize*adjBox+aMol);
              anAdjAdjTars[aBox].push_back(aTar);
            }
        }
      adjMols.resize(0);
    }
}

void Species::setRands(const unsigned currBox,
                       const unsigned r,
                       unsigned aSize,
                       const std::vector<unsigned>& anAdjBoxes,
                       std::vector<unsigned>& aRands,
                       RandomLib::Random& aRng)

{
#if defined(__MIC__)
  const unsigned BLOCKSIZE(16);
  for(unsigned iStart(0); iStart < anAdjBoxes.size(); iStart += BLOCKSIZE)
    {
      const unsigned iStop = std::min(iStart + BLOCKSIZE, (unsigned)anAdjBoxes.size());
      for(unsigned i(iStart); i < iStop; ++i)
        {
          const unsigned adjBox(anAdjBoxes[i]/(theBoxSize/theThreadSize));
          const unsigned aBoxID(anAdjBoxes[i]%(theBoxSize/theThreadSize));
          _mm_prefetch((const char*)(theThreads[adjBox]->getAdjMolsAddress(currBox, r, aBoxID)), _MM_HINT_T0);
        }
      for(unsigned i(iStart); i < iStop; ++i)
        {
          const unsigned adjBox(anAdjBoxes[i]/(theBoxSize/theThreadSize));
          const unsigned aBoxID(anAdjBoxes[i]%(theBoxSize/theThreadSize));
          aSize += theThreads[adjBox]->getAdjMolsSize(currBox, r, aBoxID);
        }
    }
  aRands.resize(aSize);
  for(unsigned i(0); i < aSize; ++i)
    {
      aRands[i] = aRng.IntegerC(11);
    }
#else
  for(unsigned i(0); i != anAdjBoxes.size(); ++i)
    {
      const unsigned adjBox(anAdjBoxes[i]/(theBoxSize/theThreadSize));
      const unsigned aBoxID(anAdjBoxes[i]%(theBoxSize/theThreadSize));
      aSize += theThreads[adjBox]->getAdjMolsSize(currBox, r, aBoxID);
    }
  aRands.resize(aSize);
  for(unsigned i(0); i < aSize; ++i)
    {
      aRands[i] = aRng.IntegerC(11);
    }
#endif
}

void Species::setTars(const unsigned currBox,
                      std::vector<unsigned>& aMols,
                      std::vector<unsigned>& aTars,
                      std::vector<unsigned>* anAdjMols,
                      std::vector<unsigned>* anAdjTars,
                      const std::vector<unsigned>& anAdjoins,
                      std::vector<unsigned>& aRands)
{
#if defined(__MIC__)
  const unsigned BLOCKSIZE(64);
  const unsigned aSize(aMols.size());
  aTars.resize(aMols.size());
  unsigned j(0);

  for(unsigned iStart(0); iStart < aSize; iStart += BLOCKSIZE)
    {
      const unsigned iStop = std::min(iStart + BLOCKSIZE, aSize);
      for(unsigned i(iStart); i < iStop; ++i)
        {
          const unsigned index(aMols[i]*theAdjoinSize+aRands[i]);
          _mm_prefetch((const char*)(anAdjoins.data() + index), _MM_HINT_T0);
        }
      for(unsigned i(iStart); i < iStop; ++i)
        {
          const unsigned aMol = aMols[i];
          const unsigned index(aMols[i]*theAdjoinSize+aRands[i]);
          const unsigned aTar = anAdjoins[index];
          if(aTar/theBoxMaxSize != currBox)
            {
              anAdjMols[aTar/theBoxMaxSize].push_back(aMol);
              anAdjTars[aTar/theBoxMaxSize].push_back(aTar);
            }
          else
            {
              aMols[j] = aMol;
              aTars[j] = aTar;
              ++j;
            }
        }
    }
  aMols.resize(j);
  aTars.resize(j);
#else
  const unsigned aSize(aMols.size());
  aTars.resize(aSize);
  for(unsigned i(0); i < aSize; ++i)
    {
      aTars[i] = anAdjoins[aMols[i]*theAdjoinSize+aRands[i]];
    }
  for(unsigned i(0); i < aMols.size(); ++i)
    {
      if(aTars[i]/theBoxMaxSize != currBox)
        {
          anAdjMols[aTars[i]/theBoxMaxSize].push_back(aMols[i]);
          anAdjTars[aTars[i]/theBoxMaxSize].push_back(aTars[i]);
          aMols[i] = aMols.back();
          aMols.pop_back();
          aTars[i] = aTars.back();
          aTars.pop_back();
          --i;
        }
    }
#endif
}

void Species::walk(const unsigned anID, unsigned r, unsigned w,
           RandomLib::Random& aRng,
           std::vector<unsigned>& aMols,
           std::vector<unsigned>& aTars,
           std::vector<unsigned>* anAdjMols,
           std::vector<unsigned>* anAdjTars,
           std::vector<std::vector<std::vector<unsigned> > >& anAdjAdjMols,
           std::vector<std::vector<std::vector<unsigned> > >& anAdjAdjTars,
           unsigned* aBorderCounts,
           std::vector<unsigned>* aBorderMols,
           std::vector<unsigned>* aBorderTars,
           std::vector<std::vector<unsigned> >& aRepeatAdjMols,
           std::vector<std::vector<unsigned> >& aRepeatAdjTars,
           const std::vector<unsigned>& anAdjoins,
           std::vector<unsigned short>& anIDs,
           const std::vector<unsigned>& anAdjBoxes,
           const std::vector<unsigned>& anAdjAdjBoxes,
           std::vector<unsigned>& aRands)
{
  unsigned aTotalBoxSize = theStepper->getBoxSize();
  unsigned aBorderOffset = aTotalBoxSize * w;
  updateBoxMols(anID, r, aMols, aTars, anAdjBoxes);
  walkMols(aMols, aTars, anIDs);
  updateAdjMols(anID, r, aRepeatAdjMols, aRepeatAdjTars, anAdjBoxes);
  updateAdjAdjMols(anID, r, anAdjAdjBoxes);
  walkAdjMols(anID, r, aMols, anIDs, anAdjBoxes);
  setRands(anID, r, aMols.size(), anAdjBoxes, aRands, aRng);
  setTars(anID, aMols, aTars, anAdjMols + aTotalBoxSize * w, anAdjTars + aTotalBoxSize * w, anAdjoins, aRands);
  setAdjTars(anID, r, aBorderCounts + aBorderOffset, aBorderMols + aBorderOffset, aBorderTars + aBorderOffset, anAdjAdjMols[w], anAdjAdjTars[w], aRepeatAdjMols, aRepeatAdjTars, anAdjBoxes, aMols.size(), aRands);
}

void Species::updateMols()
{
  if(isDiffusiveVacant || isReactiveVacant)
    {
      updateVacantMols();
    }
  else if(isTag)
    {
      updateTagMols();
    }
  if(!theID && theStepper->getCurrentTime() > 0)
    {
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          const unsigned aThreadID(i/(theBoxSize/theThreadSize));
          const unsigned aBoxID(i%(theBoxSize/theThreadSize));
          theThreads[aThreadID]->updateMols(theMols[i], aBoxID);
        }
    }
}

}
