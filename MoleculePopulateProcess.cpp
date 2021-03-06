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

#include <algorithm>
#include <gsl/gsl_randist.h>
#include <boost/lexical_cast.hpp>
#include <libecs/Model.hpp>
#include <libecs/System.hpp>
#include <libecs/Stepper.hpp>
#include <libecs/Process.hpp>
#include <libecs/VariableReference.hpp>
#include "MoleculePopulateProcess.hpp"
#include "Species.hpp"
#include "SpatiocyteProcess.hpp"

LIBECS_DM_INIT(MoleculePopulateProcess, Process); 


void MoleculePopulateProcess::initialize()
{
  if(isInitialized)
    {
      return;
    }
  SpatiocyteProcess::initialize();
  isPriorityQueued = true;
  for(VariableReferenceVector::const_iterator
	i(theVariableReferenceVector.begin());
        i != theVariableReferenceVector.end(); ++i)
    {
      Variable* aVariable((*i).getVariable());
      if(aVariable->getName() == "HD")
        {
          THROW_EXCEPTION(ValueError, getPropertyInterface().getClassName() +
            " [" + getFullID().asString() + "]: " +  
            aVariable->getFullID().asString() + " is a HD species and " +
            "therefore cannot be populated");
        }
    }
}

void MoleculePopulateProcess::initializeSecond()
{
  SpatiocyteProcess::initializeSecond();
  checkProcess();
  for(std::vector<Species*>::const_iterator i(theProcessSpecies.begin());
      i != theProcessSpecies.end(); ++i)
    {
      (*i)->setPopulateProcess(this, GaussianSigma);
    }
}

void MoleculePopulateProcess::checkProcess()
{
  std::vector<Variable*> aVariables;
  std::vector<Process*> aProcesses(theSpatiocyteStepper->getProcessVector());
  for(std::vector<Process*>::const_iterator j(aProcesses.begin());
      j!= aProcesses.end(); ++j)
    {
      MoleculePopulateProcess* aProcess(
                              dynamic_cast<MoleculePopulateProcess*>(*j));
      if(aProcess)
        {
          VariableReferenceVector aVariableReferences(
                                  (*j)->getVariableReferenceVector());
          for(VariableReferenceVector::const_iterator 
              k(aVariableReferences.begin());
              k != aVariableReferences.end(); ++k)
            {
              aVariables.push_back((*k).getVariable());
            }	
        }     
    }
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i !=theSpecies.end(); ++i)
    {
      if((*i)->getVariable() && (*i)->getVariable()->getValue() != 0 &&
         std::find(aVariables.begin(), aVariables.end(),
                   (*i)->getVariable()) == aVariables.end())
        {
          THROW_EXCEPTION(ValueError, getPropertyInterface().getClassName() +
            " [" + (*i)->getVariable()->getFullID().asString() + "]: " +  
            "has a non-zero value but is not populated.");
        }
    }
}

void MoleculePopulateProcess::fire()
{
  for(std::vector<Species*>::const_iterator i(theProcessSpecies.begin());
      i != theProcessSpecies.end(); ++i)
    {
      (*i)->removeMols();
      populateUniformSparse(*i);
    }
  theInterval = ResetTime;
  theTime += theInterval; 
  thePriorityQueue->move(theQueueID);
}

void MoleculePopulateProcess::populateGaussian(Species* aSpecies)
{
}

void MoleculePopulateProcess::populateUniformOnMultiscale(Species* aSpecies)
{
  /*
  std::cout << "    Populating uniformly on multiscale vacant:" <<
    getIDString(aSpecies) << " current size:" << aSpecies->size() <<
    ", populate size:" << aSpecies->getPopulateMolSize() << std::endl;
  if(!aSpecies->getIsPopulated())
    {
      if(UniformRadiusX == 1 && UniformRadiusY == 1 && UniformRadiusZ == 1 &&
         !OriginX && !OriginY && !OriginZ)
        {
          Species* aVacantSpecies(aSpecies->getVacantSpecies());
          aVacantSpecies->updateMols();
          unsigned int aSize(aSpecies->getPopulateMolSize());
          unsigned int aVacantSize(aVacantSpecies->getPopulatableSize());
          if(aVacantSize < aSize)
            {
              THROW_EXCEPTION(ValueError, String(
                              getPropertyInterface().getClassName()) +
                              "[" + getFullID().asString() + "]: There are " +
                              int2str(aSize) + " " + getIDString(aSpecies) +
                              " molecules that must be uniformly populated," +
                              "\nbut there are only " +
                              int2str(aVacantSize) + 
                              " multiscale vacant voxels of " +
                              getIDString(aSpecies->getVacantSpecies()) +
                              " that can be populated on.");
            }
          for(unsigned int i(0); i != aSize; ++i)
            {
              unsigned aMol;
              //After a molecule is added, the diffusive vacant species
              //molecule list is not updated, so we need to check if
              //it really is still vacant:
              aMol = aVacantSpecies->getRandomPopulatableMol();
              aSpecies->addMol(aMol);
            }
        }
      else
        {
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
    */
}

void MoleculePopulateProcess::populateUniformOnDiffusiveVacant(Species*
                                                               aSpecies)
{
  /*
  std::cout << "    Populating uniformly on diffusive vacant:" <<
    getIDString(aSpecies) << " current size:" << aSpecies->size() <<
    ", populate size:" << aSpecies->getPopulateMolSize() << std::endl;
  if(!aSpecies->getIsPopulated())
    {
      if(UniformRadiusX == 1 && UniformRadiusY == 1 && UniformRadiusZ == 1 &&
         !OriginX && !OriginY && !OriginZ)
        {
          Species* aVacantSpecies(aSpecies->getVacantSpecies());
          aVacantSpecies->updateMols();
          unsigned int aSize(aSpecies->getPopulateMolSize());
          unsigned int aVacantSize(aVacantSpecies->getPopulatableSize());
          if(aVacantSize < aSize)
            {
              THROW_EXCEPTION(ValueError, String(
                              getPropertyInterface().getClassName()) +
                              "[" + getFullID().asString() + "]: There are " +
                              int2str(aSize) + " " + getIDString(aSpecies) +
                              " molecules that must be uniformly populated," +
                              "\nbut there are only " +
                              int2str(aVacantSize) + 
                              " diffuse vacant voxels of " +
                              getIDString(aSpecies->getVacantSpecies()) +
                              " that can be populated on.");
            }
          for(unsigned int i(0); i != aSize; ++i)
            {
              unsigned aMol;
              //After a molecule is added, the diffusive vacant species
              //molecule list is not updated, so we need to check if
              //it really is still vacant:
              aMol = aVacantSpecies->getRandomPopulatableMol();
              aSpecies->addMol(aMol);
            }
        }
      else
        {
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
    */
}

void MoleculePopulateProcess::populateUniformDense(Species* aSpecies,
                                              unsigned int* aList, 
                                              unsigned int* aCount)
{
  /*
  std::cout << "    Populating densely:" <<
    getIDString(aSpecies) << " current size:" << aSpecies->size() <<
    ", populate size:" << aSpecies->getPopulateMolSize() << std::endl;
  Species* aVacantSpecies(aSpecies->getVacantSpecies());
  if(!aSpecies->getIsPopulated())
    {
      if(UniformRadiusX == 1 && UniformRadiusY == 1 && UniformRadiusZ == 1 &&
         !OriginX && !OriginY && !OriginZ)
        {
          unsigned int aSize(aSpecies->getPopulateMolSize());
          for(unsigned int j(0); j != aSize; ++j)
            {
              unsigned int aMol;
              do
                {
                  aMol = aVacantSpecies->getMol(aList[(*aCount)++]); 
                }
              while((*theIDs)[aMol] != aVacantSpecies->getID());
              aSpecies->addMol(aMol);
            }
        }
      else
        {
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
    */
}

void MoleculePopulateProcess::populateUniformSparse(Species* aSpecies)
{
  std::cout << "    Populating sparsely:" <<
    getIDString(aSpecies) << " current size:" << aSpecies->size() <<
    ", populate size:" << aSpecies->getTotalPopulateMolSize() << std::endl;
  Species* aVacantSpecies(aSpecies->getVacantSpecies());
  if(!aSpecies->getIsPopulated())
    {
      if(UniformRadiusX == 1 && UniformRadiusY == 1 && UniformRadiusZ == 1 &&
         !OriginX && !OriginY && !OriginZ)
        {
          for(unsigned i(0); i != theIDs->size(); ++i)
          //for(unsigned i(0); i != 1;  ++i)
            {
              unsigned int aSize(aSpecies->getPopulateMolSize(i));
              int availableVoxelSize(aVacantSpecies->size(i));
              for(unsigned int j(0); j != aSize; ++j)
                {
                  unsigned int aMol;
                  do
                    {
                      aMol = aVacantSpecies->getMol(i, gsl_rng_uniform_int(
                                getStepper()->getRng(), availableVoxelSize));
                    }
                  while((*theIDs)[i][aMol] != aVacantSpecies->getID());
                  aSpecies->addMol(i, aMol);
                }
            }
        }
      else
        { 
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
  aSpecies->updateMols();
}

void MoleculePopulateProcess::populateUniformRanged(Species* aSpecies)
{
  /*
  std::cout << "    Populating uniformly ranged:" <<
    getIDString(aSpecies) << " current size:" << aSpecies->size() <<
    ", populate size:" << aSpecies->getPopulateMolSize() << std::endl;
  Comp* aComp(aSpecies->getComp());
  Species* aVacantSpecies(aSpecies->getVacantSpecies());
  double deltaX(0);
  double deltaY(0);
  double deltaZ(0);
  // Increase the compartment dimensions by delta if it is a surface 
  // compartment:
  if(aComp->dimension == 2)
    {
      deltaX = theSpatiocyteStepper->getNormalizedVoxelRadius()*6/
        aComp->lengthX;
      deltaY = theSpatiocyteStepper->getNormalizedVoxelRadius()*6/
        aComp->lengthY;
      deltaZ = theSpatiocyteStepper->getNormalizedVoxelRadius()*6/
        aComp->lengthZ;
    }
  double maxX(std::min(1.0, OriginX+UniformRadiusX));
  double minX(std::max(-1.0, OriginX-UniformRadiusX));
  double maxY(std::min(1.0, OriginY+UniformRadiusY));
  double minY(std::max(-1.0, OriginY-UniformRadiusY));
  double maxZ(std::min(1.0, OriginZ+UniformRadiusZ));
  double minZ(std::max(-1.0, OriginZ-UniformRadiusZ)); 
  maxX = aComp->centerPoint.x + maxX*aComp->lengthX/2*(1+deltaX);
  minX = aComp->centerPoint.x + minX*aComp->lengthX/2*(1+deltaX);
  maxY = aComp->centerPoint.y + maxY*aComp->lengthY/2*(1+deltaY);
  minY = aComp->centerPoint.y + minY*aComp->lengthY/2*(1+deltaY);
  maxZ = aComp->centerPoint.z + maxZ*aComp->lengthZ/2*(1+deltaZ);
  minZ = aComp->centerPoint.z + minZ*aComp->lengthZ/2*(1+deltaZ);
  std::vector<unsigned int> aMols;
  unsigned int aVacantSize(aVacantSpecies->getPopulatableSize());
  for(unsigned int i(0); i != aVacantSize; ++i)
    {
      unsigned int aMol(aVacantSpecies->getPopulatableMol(i));
      Point& aPoint((*theInfo)[aMol].point);
      if((*theIDs)[aMol] == aSpecies->getVacantID() &&
         aPoint.x < maxX && aPoint.x > minX &&
         aPoint.y < maxY && aPoint.y > minY &&
         aPoint.z < maxZ && aPoint.z > minZ)
        {
          aMols.push_back(aMol);
        }
    }
  unsigned int aSize(aSpecies->getPopulateMolSize());
  if(aMols.size() < aSize)
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + "]: There are " +
                      int2str(aSize) + " " + getIDString(aSpecies) +
                      " molecules that must be uniformly populated in a " +
                      "given range,\n but there are only " +
                      int2str(aMols.size()) + " vacant voxels of " +
                      getIDString(aSpecies->getVacantSpecies()) +
                      " that can be populated.");
    }
  std::random_shuffle(aMols.begin(), aMols.end());
  for(unsigned int i(0); i != aMols.size(); ++i)
    {
      unsigned aMol(aMols[i]);
      if(aSpecies->size() < aSize && aSpecies->isPopulatable(aMol))
        {
          aSpecies->addMol(aMol);
        }
    }
    */
}

