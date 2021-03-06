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

#include "CompartmentProcess.hpp"
#include "Vector.hpp"

LIBECS_DM_INIT(CompartmentProcess, Process); 

unsigned CompartmentProcess::getLatticeResizeMol(unsigned aStartMol)
{
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  aComp->interfaceID = theInterfaceSpecies->getID();
  *theComp = *aComp;
  theVacantSpecies->resetFixedAdjoins();
  theVacantSpecies->setMolRadius(DiffuseRadius);
  if(theLipidSpecies)
    {
      theLipidSpecies->resetFixedAdjoins();
      theLipidSpecies->setMolRadius(LipidRadius);
    }
  //The compartment center point (origin):
  Origin = aComp->centerPoint;
  Origin.x += OriginX*aComp->lengthX/2;
  Origin.y += OriginY*aComp->lengthY/2;
  Origin.z += OriginZ*aComp->lengthZ/2;
  setCompartmentDimension();
  theComp->dimension = dimension;
  setLipidCompSpeciesProperties();
  setVacantCompSpeciesProperties();
  subStartMol = aStartMol;
  lipStartMol = aStartMol+Filaments*Subunits;
  endMol = lipStartMol+LipidRows*LipidCols;
  return endMol-aStartMol;
}

void CompartmentProcess::setVacantCompSpeciesProperties()
{
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      theVacantCompSpecies[i]->setDimension(dimension);
      theVacantCompSpecies[i]->setMolRadius(SubunitRadius);
      theVacantCompSpecies[i]->setDiffuseRadius(DiffuseRadius);
      if(theLipidSpecies)
        {
          theVacantCompSpecies[i]->setMultiscaleVacantSpecies(theLipidSpecies);
        }
    }
}

void CompartmentProcess::setVacantCompMultiscaleProperties()
{
  for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
    {
      Species* aLipid(theLipidCompSpecies[i]);
      for(unsigned j(0); j != theVacantCompSpecies.size(); ++j)
        {
          Species* aVacant(theVacantCompSpecies[j]);
          if(aLipid->getVacantSpecies() == aVacant)
            {
              int aCoeff(getCoefficient(aLipid));
              Species* aLipidPair(coefficient2species(sqrt(aCoeff*aCoeff)));
              aVacant->setMultiscaleBindUnbindIDs(aLipid->getID(),
                                                  aLipidPair->getID());
            }
        }
    }
}

int CompartmentProcess::getCoefficient(Species* aSpecies)
{
  for(VariableReferenceVector::iterator i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      if(aSpecies->getVariable() == (*i).getVariable()) 
        {
          return (*i).getCoefficient();
        }
    }
  return 0;
}

Species* CompartmentProcess::coefficient2species(int aCoeff)
{
  for(VariableReferenceVector::iterator i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      if((*i).getCoefficient() == aCoeff)
        {
          return theSpatiocyteStepper->variable2species((*i).getVariable());
        }
    }
  return NULL;
}

void CompartmentProcess::setLipidCompSpeciesProperties()
{
  for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
    {
      theLipidCompSpecies[i]->setDimension(dimension);
      theLipidCompSpecies[i]->setMolRadius(LipidRadius);
    }
}

void CompartmentProcess::updateResizedLattice()
{
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      //TODO: replace subStartMol with theVacantSpecies->getMol(0) 
      //TODO: replace lipStartMol with theLipidSpecies->getMol(0) 
      theVacantCompSpecies[i]->setVacStartMol(subStartMol);
      theVacantCompSpecies[i]->setLipStartMol(lipStartMol);
    }
}

void CompartmentProcess::setCompartmentDimension()
{
  if(Length)
    {
      Subunits = (unsigned)rint(Length/(DiffuseRadius*2));
    }
  if(Width)
    {
      Filaments = (unsigned)rint((Width-2*DiffuseRadius)/
                                 (DiffuseRadius*sqrt(3)))+1;
    }
  if(Periodic && Filaments%2 != 0)
    {
      ++Filaments;
    }
  Width = 2*DiffuseRadius+(Filaments-1)*DiffuseRadius*sqrt(3); 
  Height = 2*DiffuseRadius;
  if(Filaments == 1)
    {
      dimension = 1;
      Length = Subunits*DiffuseRadius*2;
    }
  else
    {
      dimension = 2;
      //Add DiffuseRadius for the protrusion from hexagonal arrangement:
      Length = Subunits*DiffuseRadius*2+DiffuseRadius;
    }
  //Normalized compartment lengths in terms of lattice voxel radius:
  nLength = Length/(VoxelRadius*2);
  nWidth = Width/(VoxelRadius*2);
  nHeight = Height/(VoxelRadius*2);
  theComp->lengthX = nHeight;
  theComp->lengthY = nLength;
  theComp->lengthZ = nWidth;
  if(theLipidSpecies)
    {
      LipidCols = (unsigned)rint(Length/(LipidRadius*2));
      LipidRows = (unsigned)rint((Width-2*LipidRadius)/(LipidRadius*sqrt(3)))+1;
    }
  //Actual surface area = Width*Length
}

void CompartmentProcess::initializeThird()
{
  if(!isCompartmentalized)
    {
      //setVacantCompMultiscaleProperties must be in 
      //initializeThird since it requires vacant species properties
      //set by DiffusionProcess in initializeSecond:
      setVacantCompMultiscaleProperties();

      initializeVectors();
      initializeFilaments(subunitStart, Filaments, Subunits, nDiffuseRadius,
                          theVacantSpecies, subStartMol);
      elongateFilaments(theVacantSpecies, subStartMol, Filaments, Subunits,
                        nDiffuseRadius);
      connectFilaments(subStartMol, Filaments, Subunits);
      interfaceSubunits();
      initializeFilaments(lipidStart, LipidRows, LipidCols, nLipidRadius,
                          theLipidSpecies, lipStartMol);
      elongateFilaments(theLipidSpecies, lipStartMol, LipidRows,
                        LipidCols, nLipidRadius);
      connectFilaments(lipStartMol, LipidRows, LipidCols);
      setDiffuseSize(lipStartMol, endMol);
      setSpeciesIntersectLipids();
      isCompartmentalized = true;
    }
  theVacantSpecies->setIsPopulated();
  theLipidSpecies->setIsPopulated();
}

void CompartmentProcess::setSpeciesIntersectLipids()
{
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      theVacantCompSpecies[i]->setIntersectLipids(theLipidSpecies);
    }
}

Point CompartmentProcess::getStartVoxelPoint()
{
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  Species* surface(aComp->surfaceSub->vacantSpecies);
  Point nearest;
  Point origin;
  origin.x = 0;
  origin.y = 0;
  origin.z = 0;
  double dist;
  if(surface->size())
    {
      nearest = surface->getPoint(0);
      dist = getDistance(nearest, origin);
    }
  for(unsigned i(1); i < surface->size(); ++i)
    {
      Point& aPoint(surface->getPoint(i));
      double aDist(getDistance(aPoint, origin));
      if(aDist < dist)
        {
          dist = aDist;
          nearest = aPoint;
        }
    }
  return nearest;
}

void CompartmentProcess::initializeVectors()
{
  subunitStart = getStartVoxelPoint();
  lengthStart = subunitStart;
  lengthStart.z -= 2*nDiffuseRadius;
  lengthStart.y -= nDiffuseRadius;

  lengthVector.x = 0;
  lengthVector.y = 0;
  lengthVector.z = 1;
  lengthEnd = disp(lengthStart, lengthVector, nLength);

  widthVector.x = 0;
  widthVector.y = 1;
  widthVector.x = 0;
  widthEnd = disp(lengthEnd, widthVector, nWidth);

  heightVector.x = 1;
  heightVector.y = 0;
  heightVector.z = 0;
  heightEnd = disp(widthEnd, heightVector, nHeight);

  if(theLipidSpecies)
    {
      lipidStart = lengthStart;
      disp_(lipidStart, lengthVector, nLipidRadius);
      disp_(lipidStart, widthVector, nLipidRadius);
    }

  Point center(lengthStart);
  disp_(center, lengthVector, nLength/2);
  disp_(center, widthVector, nWidth/2);
  theComp->centerPoint = center;

  //Set up surface vectors:
  surfaceNormal = cross(lengthVector, widthVector);
  surfaceNormal = norm(surfaceNormal);
  surfaceDisplace = dot(surfaceNormal, widthEnd);
  lengthDisplace = dot(lengthVector, lengthStart);
  lengthDisplaceOpp = dot(lengthVector, lengthEnd);
  widthDisplace = dot(widthVector, widthEnd);
  widthDisplaceOpp = dot(widthVector, lengthEnd);
}

void CompartmentProcess::rotate(Point& V)
{
  theSpatiocyteStepper->rotateX(RotateX, &V, -1);
  theSpatiocyteStepper->rotateY(RotateY, &V, -1);
  theSpatiocyteStepper->rotateZ(RotateZ, &V, -1);
}

void CompartmentProcess::initializeFilaments(Point& aStartPoint, unsigned aRows,
                                             unsigned aCols, double aRadius,
                                             Species* aVacant,
                                             unsigned aStartMol)
{
  //The first comp voxel must have the aStartMol:
  addCompVoxel(0, 0, aStartPoint, aVacant, aStartMol, aCols);
  for(unsigned i(1); i != aRows; ++i)
    {
      Point U(aStartPoint);
      disp_(U, widthVector, i*aRadius*sqrt(3)); 
      if(i%2 == 1)
        {
          disp_(U, lengthVector, -aRadius); 
        }
      addCompVoxel(i, 0, U, aVacant, aStartMol, aCols);
    }
}

void CompartmentProcess::addCompVoxel(unsigned rowIndex, 
                                      unsigned colIndex,
                                      Point& aPoint,
                                      Species* aVacant,
                                      unsigned aStartMol,
                                      unsigned aCols)
{
  unsigned aMol(aStartMol+rowIndex*aCols+colIndex);
  VoxelInfo& anInfo((*theInfo)[aMol]);
  anInfo.point = aPoint;
  anInfo.adjoinSize = 0;
  aVacant->addCompVoxel(aMol);
}

void CompartmentProcess::elongateFilaments(Species* aVacant,
                                           unsigned aStartMol,
                                           unsigned aRows,
                                           unsigned aCols,
                                           double aRadius)
{
  for(unsigned i(0); i != aRows; ++i)
    {
      VoxelInfo& anInfo((*theInfo)[aStartMol+i*aCols]);
      Point A(anInfo.point);
      for(unsigned j(1); j != aCols; ++j)
        {
          disp_(A, lengthVector, aRadius*2);
          addCompVoxel(i, j, A, aVacant, aStartMol, aCols);
        }
    }
}

/*
 [NW] [  N  ] [NE]
     [subunit] 
 [SW] [  S  ] [SE]
 */

void CompartmentProcess::connectFilaments(unsigned aStartMol,
                                          unsigned aRows, unsigned aCols)
{
  for(unsigned i(0); i != aCols; ++i)
    {
      for(unsigned j(0); j != aRows; ++j)
        { 
          if(i > 0)
            { 
              //N-S
              unsigned a(aStartMol+j*aCols+(i-1));
              unsigned b(aStartMol+j*aCols+i);
              connectSubunit(a, b);
            }
          else if(Periodic)
            {
              //periodic N-S
              unsigned a(aStartMol+j*aCols); 
              unsigned b(aStartMol+j*aCols+aCols-1);
              connectSubunit(a, b);
              if(j%2 == 1)
                {
                  if(j+1 < aRows)
                    {
                      //periodic NE_SW 
                      b = aStartMol+(j+1)*aCols+aCols-1; 
                      connectSubunit(a, b);
                    }
                  else if(j == aRows-1)
                    {
                      //periodic NE_SW 
                      b = aStartMol+aCols-1; 
                      connectSubunit(a, b);
                    }
                  //periodic SE_NW 
                  b = aStartMol+(j-1)*aCols+aCols-1; 
                  connectSubunit(a, b);
                }
            }
          if(j > 0)
            {
              if(j%2 == 1)
                {
                  //NE_SW
                  unsigned a(aStartMol+j*aCols+i);
                  unsigned b(aStartMol+(j-1)*aCols+i); 
                  connectSubunit(a, b);
                  if(i > 0)
                    {
                      //SE_NW
                      b = aStartMol+(j-1)*aCols+(i-1); 
                      connectSubunit(a, b);
                    }
                }
              else
                {
                  //SE_NW
                  unsigned a(aStartMol+j*aCols+i);
                  unsigned b(aStartMol+(j-1)*aCols+i); 
                  connectSubunit(a, b);
                  if(i+1 < aCols)
                    {
                      //NE_SW
                      b = aStartMol+(j-1)*aCols+(i+1); 
                      connectSubunit(a, b);
                    }
                }
            }
        }
      if(Periodic && aRows > 1)
        { 
          //periodic SE_NW
          unsigned a(aStartMol+i); //row 0
          unsigned b(aStartMol+(aRows-1)*aCols+i); 
          connectSubunit(a, b);
          if(i+1 < aCols)
            {
              //periodic NE_SW
              b = aStartMol+(aRows-1)*aCols+(i+1); 
              connectSubunit(a, b);
            }
        }
    }
}

void CompartmentProcess::setDiffuseSize(unsigned start, unsigned end)
{
  for(unsigned i(start); i != end; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
      subunit.diffuseSize = (*theInfo)[i].adjoinSize;
    }
}

void CompartmentProcess::connectSubunit(unsigned a, unsigned b)
{
  addAdjoin(a, b);
  addAdjoin(b, a);
}

void CompartmentProcess::addAdjoin(unsigned sub, unsigned mol)
{
  VoxelInfo& anInfo((*theInfo)[sub]);
  Voxel& aVoxel((*theLattice)[sub]);
  unsigned* temp(new unsigned[anInfo.adjoinSize+1]);
  for(unsigned int i(0); i != anInfo.adjoinSize; ++i)
    {
      //Avoid duplicated adjoins:
      if(aVoxel.adjoins[i] == mol)
        {
          delete[] temp;
          return;
        }
      temp[i] = aVoxel.adjoins[i];
    }
  delete[] aVoxel.adjoins;
  temp[anInfo.adjoinSize++] = mol;
  aVoxel.adjoins = temp;
}

void CompartmentProcess::interfaceSubunits()
{
  enlistInterfaceVoxels();
  enlistNonIntersectInterfaceVoxels();
  setDiffuseSize(subStartMol, lipStartMol);
  enlistSubunitInterfaceAdjoins();
  theVacantSpecies->setIsPopulated();
  theInterfaceSpecies->setIsPopulated();
}

void CompartmentProcess::enlistInterfaceVoxels()
{
  subunitInterfaces.resize(Filaments*Subunits);
  for(unsigned i(subStartMol); i != lipStartMol; ++i)
    {
      Point bottomLeft((*theInfo)[i].point);
      Point topRight((*theInfo)[i].point);
      bottomLeft.x -= nSubunitRadius+theSpatiocyteStepper->getColLength();
      bottomLeft.y -= nSubunitRadius+theSpatiocyteStepper->getLayerLength();
      bottomLeft.z -= nSubunitRadius+theSpatiocyteStepper->getRowLength();
      topRight.x += nSubunitRadius+theSpatiocyteStepper->getColLength();
      topRight.y += nSubunitRadius+theSpatiocyteStepper->getLayerLength();
      topRight.z += nSubunitRadius+theSpatiocyteStepper->getRowLength();
      unsigned blRow(0);
      unsigned blLayer(0);
      unsigned blCol(0);
      theSpatiocyteStepper->point2global(bottomLeft, blRow, blLayer, blCol);
      unsigned trRow(0);
      unsigned trLayer(0);
      unsigned trCol(0);
      theSpatiocyteStepper->point2global(topRight, trRow, trLayer, trCol);
      for(unsigned j(blRow); j <= trRow; ++j)
        {
          for(unsigned k(blLayer); k <= trLayer; ++k)
            {
              for(unsigned l(blCol); l <= trCol; ++l)
                {
                  unsigned m(theSpatiocyteStepper->global2coord(j, k, l));
                  addInterfaceVoxel(i, m);
                }
            }
        }
    }
}

void CompartmentProcess::addInterfaceVoxel(unsigned subunitMol,
                                           unsigned voxelMol)
{ 
  Point& subunitPoint((*theInfo)[subunitMol].point);
  Point voxelPoint(theSpatiocyteStepper->coord2point(voxelMol));
  double dist(getDistance(subunitPoint, voxelPoint));
  //Should use SubunitRadius instead of DiffuseRadius since it is the
  //actual size of the subunit:
  if(dist <= (nSubunitRadius+nVoxelRadius)*1.0001) 
    {
      Voxel& voxel((*theLattice)[voxelMol]);
      //theSpecies[6]->addMol(&voxel);
      //Insert voxel in the list of interface voxels if was not already:
      //if(voxel.id == theComp->vacantSpecies->getID())
      if(voxel.id != theInterfaceSpecies->getID())
        {
          //theSpecies[5]->addMol(&voxel);
          theInterfaceSpecies->addMol(voxelMol);
        }
      //each subunit list has unique (no duplicates) interface voxels:
      subunitInterfaces[subunitMol-subStartMol].push_back(voxelMol);
    }
}

void CompartmentProcess::enlistSubunitInterfaceAdjoins()
{
  for(unsigned i(0); i != subunitInterfaces.size(); ++i)
    {
      for(unsigned j(0); j != subunitInterfaces[i].size(); ++j)
        {
          addAdjoin(subunitInterfaces[i][j], i+subStartMol);
          Voxel& interface((*theLattice)[subunitInterfaces[i][j]]);
          for(unsigned k(0); k != interface.diffuseSize; ++k)
            {
              unsigned mol(interface.adjoins[k]);
              Voxel& adjoin((*theLattice)[mol]);
              if(adjoin.id != theInterfaceSpecies->getID())
                {
                  addAdjoin(i+subStartMol, mol);
                }
            }
        }
    }
}

void CompartmentProcess::enlistNonIntersectInterfaceVoxels()
{
  for(unsigned i(0); i != theInterfaceSpecies->size(); ++i)
    {
      unsigned voxelMol(theInterfaceSpecies->getMol(i));
      Voxel& anInterface((*theLattice)[voxelMol]);
      for(unsigned j(0); j != theAdjoinSize; ++j)
        {
          unsigned adjoin(anInterface.adjoins[j]);
          if((*theLattice)[adjoin].id != theInterfaceSpecies->getID())
            {
              Point aPoint(theSpatiocyteStepper->coord2point(adjoin));
              if(isInside(aPoint))
                {
                  addNonIntersectInterfaceVoxel(adjoin, aPoint);
                }
            }
        }
    }
}

bool CompartmentProcess::isInside(Point& aPoint)
{
  double dist(point2planeDist(aPoint, lengthVector, lengthDisplace));
  if(dist >= 0)
    {
      dist = point2planeDist(aPoint, lengthVector, lengthDisplaceOpp);
      if(dist <= 0)
        {
          dist = point2planeDist(aPoint, widthVector, widthDisplaceOpp);
          if(dist <=0)
            {
              dist = point2planeDist(aPoint, widthVector, widthDisplace);
              if(dist >= 0)
                {
                  return true;
                }
            }
        }
    }
  return false;
}

void CompartmentProcess::addNonIntersectInterfaceVoxel(unsigned aMol,
                                                       Point& aPoint)
{
  Voxel& aVoxel((*theLattice)[aMol]);
  double distA(point2planeDist(aPoint, surfaceNormal, surfaceDisplace));
  for(unsigned i(0); i != theAdjoinSize; ++i)
    {
      Voxel& adjoin((*theLattice)[aVoxel.adjoins[i]]);
      if(adjoin.id != theInterfaceSpecies->getID())
        {
          Point& pointB((*theInfo)[aVoxel.adjoins[i]].point);
          double distB(point2planeDist(pointB, surfaceNormal, surfaceDisplace));
          //if not on the same side of the plane:
          if((distA < 0) != (distB < 0))
            {
              if(abs(distA) < abs(distB))
                { 
                  //theSpecies[6]->addMol(&aVoxel);
                  theInterfaceSpecies->addMol(aMol);
                  return;
                }
              else
                {
                  //theSpecies[6]->addMol(&adjoin);
                  theInterfaceSpecies->addMol(aMol);
                }
            }
        }
    }
}

void CompartmentProcess::printParameters()
{
  std::cout << getPropertyInterface().getClassName() << "[" <<
    getFullID().asString() << "]" << std::endl;
  std::cout << "  width:" << Width << " length:" << Length <<
    " area:" << Width*Length << std::endl;
  if(theLipidSpecies)
    {
      std::cout << "  " << getIDString(theLipidSpecies) << 
        " number:" << theLipidSpecies->size() << std::endl;
      for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
        {
          std::cout << "    " << getIDString(theLipidCompSpecies[i]) <<
            " number:" << theLipidCompSpecies[i]->size() << std::endl;
        }
    } 
  std::cout << "  " << getIDString(theVacantSpecies) << 
    " number:" << theVacantSpecies->size() << std::endl;
      for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
        {
          std::cout << "    " << getIDString(theVacantCompSpecies[i]) <<
            " number:" << theVacantCompSpecies[i]->size() << std::endl;
        }
}
