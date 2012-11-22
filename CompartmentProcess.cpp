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

unsigned CompartmentProcess::getLatticeResizeCoord(unsigned aStartCoord)
{
  theComp = theSpatiocyteStepper->system2Comp(getSuperSystem());
  theVacantSpecies->resetFixedAdjoins();
  theVacantSpecies->setRadius(SubunitRadius);
  if(LipidRadius)
    {
      theLipidSpecies->setRadius(LipidRadius);
    }
  //The compartment center point (origin):
  Origin = theComp->centerPoint;
  Origin.x += OriginX*theComp->lengthX/2;
  Origin.y += OriginY*theComp->lengthY/2;
  Origin.z += OriginZ*theComp->lengthZ/2;
  setCompartmentDimension();
  for(unsigned i(0); i != theCompartmentSpecies.size(); ++i)
    {
      theCompartmentSpecies[i]->setIsOffLattice();
      theCompartmentSpecies[i]->setDimension(dimension);
      theCompartmentSpecies[i]->setVacantSpecies(theVacantSpecies);
      theCompartmentSpecies[i]->setRadius(SubunitRadius);
    }
  subStartCoord = aStartCoord;
  lipStartCoord = aStartCoord+Filaments*Subunits;
  endCoord = lipStartCoord+LipidRows*LipidCols;
  return endCoord-aStartCoord;
}

void CompartmentProcess::setCompartmentDimension()
{
  if(Length)
    {
      Subunits = (unsigned)rint(Length/(SubunitRadius*2));
    }
  if(Width)
    {
      Filaments = (unsigned)rint((Width-2*SubunitRadius)/
                                 (SubunitRadius*sqrt(3)))+1;
    }
  Width = 2*SubunitRadius+(Filaments-1)*SubunitRadius*sqrt(3); 
  Height = 2*SubunitRadius;
  if(Filaments == 1)
    {
      dimension = 1;
      Length = Subunits*SubunitRadius*2;
    }
  else
    {
      dimension = 2;
      //Add SubunitRadius for the protrusion from hexagonal arrangement:
      Length = Subunits*SubunitRadius*2+SubunitRadius;
    }
  //Normalized compartment lengths in terms of lattice voxel radius:
  nLength = Length/(VoxelRadius*2);
  nWidth = Width/(VoxelRadius*2);
  nHeight = Height/(VoxelRadius*2);
  if(LipidRadius)
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
      thePoints.resize(endCoord-subStartCoord);
      initializeVectors();
      initializeFilaments(subunitStart, Filaments, Subunits, nSubunitRadius,
                          theVacantSpecies, subStartCoord);
      elongateFilaments(theVacantSpecies, subStartCoord, Filaments, Subunits,
                        nSubunitRadius);
      connectFilaments(subStartCoord, Filaments, Subunits);
      interfaceSubunits();
      initializeFilaments(lipidStart, LipidRows, LipidCols, nLipidRadius,
                          theLipidSpecies, lipStartCoord);
      elongateFilaments(theLipidSpecies, lipStartCoord, LipidRows,
                        LipidCols, nLipidRadius);
      connectFilaments(lipStartCoord, LipidRows, LipidCols);
      isCompartmentalized = true;
    }
  theVacantSpecies->setIsPopulated();
  theLipidSpecies->setIsPopulated();
}

Point CompartmentProcess::getStartVoxelPoint()
{
  Species* surface(theComp->surfaceSub->vacantSpecies);
  Point nearest;
  Point origin;
  origin.x = 0;
  origin.y = 0;
  origin.z = 0;
  double dist;
  if(surface->size())
    {
      nearest = theSpatiocyteStepper->coord2point(surface->getCoord(0));
      dist = getDistance(&nearest, &origin);
    }
  for(unsigned i(1); i < surface->size(); ++i)
    {
      Point aPoint(theSpatiocyteStepper->coord2point(surface->getCoord(i)));
      double aDist(getDistance(&aPoint, &origin));
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
  lengthStart.z -= 2*nSubunitRadius;
  lengthStart.y -= nSubunitRadius;

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

  if(LipidRadius)
    {
      lipidStart = lengthStart;
      disp_(lipidStart, lengthVector, nLipidRadius);
      disp_(lipidStart, widthVector, nLipidRadius);
    }

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
                                             unsigned aStartCoord)
{
  addCompVoxel(0, 0, aStartPoint, aVacant, aStartCoord, aCols);
  for(unsigned i(1); i != aRows; ++i)
    {
      Point U(aStartPoint);
      disp_(U, widthVector, i*aRadius*sqrt(3)); 
      if(i%2 == 1)
        {
          disp_(U, lengthVector, -aRadius); 
        }
      addCompVoxel(i, 0, U, aVacant, aStartCoord, aCols);
    }
}

void CompartmentProcess::addCompVoxel(unsigned rowIndex, 
                                      unsigned colIndex,
                                      Point& aPoint,
                                      Species* aVacant,
                                      unsigned aStartCoord,
                                      unsigned aCols)
{
  unsigned aCoord(aStartCoord+rowIndex*aCols+colIndex);
  Voxel& aVoxel((*theLattice)[aCoord]);
  aVoxel.point = &thePoints[aStartCoord-subStartCoord+rowIndex*aCols+colIndex];
  *aVoxel.point = aPoint;
  aVoxel.adjoiningCoords = new unsigned[theAdjoiningCoordSize];
  aVoxel.diffuseSize = 2;
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      aVoxel.adjoiningCoords[i] = theNullCoord;
    }
  aVacant->addCompVoxel(aCoord);
}

void CompartmentProcess::elongateFilaments(Species* aVacant,
                                           unsigned aStartCoord,
                                           unsigned aRows,
                                           unsigned aCols,
                                           double aRadius)
{
  for(unsigned i(0); i != aRows; ++i)
    {
      Voxel* startVoxel(&(*theLattice)[aStartCoord+i*aCols]);
      Point A(*startVoxel->point);
      for(unsigned j(1); j != aCols; ++j)
        {
          disp_(A, lengthVector, aRadius*2);
          addCompVoxel(i, j, A, aVacant, aStartCoord, aCols);
        }
    }
}

void CompartmentProcess::connectFilaments(unsigned aStartCoord,
                                          unsigned aRows, unsigned aCols)
{
  for(unsigned i(0); i != aCols; ++i)
    {
      for(unsigned j(0); j != aRows; ++j)
        { 
          if(i > 0)
            { 
              connectNorthSouth(aStartCoord, aCols, i, j);
            }
          else if(Periodic)
            {
              connectPeriodic(aStartCoord, aCols, j);
            }
          if(j > 0)
            {
              connectEastWest(aStartCoord, aCols, i, j);
            }
        }
      /*
      if(Filaments > 2)
        {
          connectSeamEastWest(aStartCoord, aRows, aCols, i);
        }
        */
      if(i > 0)
        {
          connectNwSw(aStartCoord, aRows, aCols, i);
        }
    }
}

void CompartmentProcess::connectPeriodic(unsigned aStartCoord, unsigned aCols,
                                         unsigned j)
{
  unsigned a(aStartCoord+j*aCols+aCols-1);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(aStartCoord+j*aCols); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[NORTH] = b;
  adjoin.adjoiningCoords[SOUTH] = a;
  aVoxel.adjoiningSize = 2;
  adjoin.adjoiningSize = 2;
}

void CompartmentProcess::connectNorthSouth(unsigned aStartCoord, unsigned aCols,
                                           unsigned i, unsigned j)
{
  unsigned a(aStartCoord+j*aCols+(i-1));
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(aStartCoord+j*aCols+i);
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[NORTH] = b;
  adjoin.adjoiningCoords[SOUTH] = a;
  aVoxel.adjoiningSize = 2;
  adjoin.adjoiningSize = 2;
}

void CompartmentProcess::connectEastWest(unsigned aStartCoord, unsigned aCols,
                                         unsigned i, unsigned j)
{
  unsigned a(aStartCoord+j*aCols+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(aStartCoord+(j-1)*aCols+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void CompartmentProcess::connectSeamEastWest(unsigned aStartCoord,
                                             unsigned aRows, unsigned aCols,
                                             unsigned i)
{
  unsigned a(aStartCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(aStartCoord+(aRows-1)*aCols+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void CompartmentProcess::connectNwSw(unsigned aStartCoord, unsigned aRows,
                                     unsigned aCols, unsigned i)
{
  unsigned a(aStartCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(aStartCoord+(aRows-1)*aCols+(i-1)); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void CompartmentProcess::interfaceSubunits()
{
  enlistInterfaceVoxels();
  enlistNonIntersectInterfaceVoxels();
  theVacantSpecies->setIsPopulated();
  theInterfaceSpecies->setIsPopulated();
}

void CompartmentProcess::enlistInterfaceVoxels()
{
  subunitInterfaces.resize(Filaments*Subunits);
  for(unsigned i(subStartCoord); i != lipStartCoord; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
      subunit.diffuseSize = subunit.adjoiningSize;
      Point center(*subunit.point);
      Point bottomLeft(*subunit.point);
      Point topRight(*subunit.point);
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

void CompartmentProcess::addInterfaceVoxel(unsigned subunitCoord,
                                           unsigned voxelCoord)
{ 
  Voxel& subunit((*theLattice)[subunitCoord]);
  Point subunitPoint(*subunit.point);
  Point voxelPoint(theSpatiocyteStepper->coord2point(voxelCoord));
  double dist(getDistance(&subunitPoint, &voxelPoint));
  if(dist <= (nSubunitRadius+nVoxelRadius)*1.0001) 
    {
      Voxel& voxel((*theLattice)[voxelCoord]);
      //theSpecies[6]->addMolecule(&voxel);
      //Insert voxel in the list of interface voxels if was not already:
      //if(voxel.id == theComp->vacantSpecies->getID())
      if(voxel.id != theInterfaceSpecies->getID())
        {
          //theSpecies[5]->addMolecule(&voxel);
          theInterfaceSpecies->addMolecule(&voxel);
        }
      subunitInterfaces[subunitCoord-subStartCoord].push_back(voxelCoord);
    }
}

void CompartmentProcess::enlistNonIntersectInterfaceVoxels()
{
  for(unsigned i(0); i != theInterfaceSpecies->size(); ++i)
    {
      unsigned voxelCoord(theInterfaceSpecies->getCoord(i));
      Voxel& anInterface((*theLattice)[voxelCoord]);
      for(unsigned j(0); j != theAdjoiningCoordSize; ++j)
        {
          Voxel& adjoin((*theLattice)[anInterface.adjoiningCoords[j]]);
          if(adjoin.id != theInterfaceSpecies->getID())
            {
              Point aPoint(theSpatiocyteStepper->coord2point(adjoin.coord));
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

void CompartmentProcess::addNonIntersectInterfaceVoxel(Voxel& aVoxel,
                                                       Point& aPoint)
{
  double distA(point2planeDist(aPoint, surfaceNormal, surfaceDisplace));
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      Voxel& adjoin((*theLattice)[aVoxel.adjoiningCoords[i]]);
      if(adjoin.id != theInterfaceSpecies->getID())
        {
          Point pointB(theSpatiocyteStepper->coord2point(adjoin.coord));
          double distB(point2planeDist(pointB, surfaceNormal, surfaceDisplace));
          //if not on the same side of the plane:
          if((distA < 0) != (distB < 0))
            {
              if(abs(distA) < abs(distB))
                { 
                  //theSpecies[6]->addMolecule(&aVoxel);
                  theInterfaceSpecies->addMolecule(&aVoxel);
                  return;
                }
              else
                {
                  //theSpecies[6]->addMolecule(&adjoin);
                  theInterfaceSpecies->addMolecule(&adjoin);
                }
            }
        }
    }
}

