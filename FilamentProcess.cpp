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

#include "FilamentProcess.hpp"

LIBECS_DM_INIT(FilamentProcess, Process); 

unsigned int FilamentProcess::getLatticeResizeCoord(unsigned int aStartCoord)
{
  theComp = theSpatiocyteStepper->system2Comp(getSuperSystem());
  theVacantSpecies->resetFixedAdjoins();
  theVacantSpecies->setRadius(VoxelRadius);
  tempID = theSpecies.size();
  C = theComp->centerPoint;
  C.x += OriginX*theComp->lengthX/2;
  C.y += OriginY*theComp->lengthY/2;
  C.z += OriginZ*theComp->lengthZ/2;
  for(unsigned int i(0); i != theFilamentSpecies.size(); ++i)
    {
      theFilamentSpecies[i]->setIsOffLattice();
      theFilamentSpecies[i]->setDimension(1);
      theFilamentSpecies[i]->setVacantSpecies(theVacantSpecies);
      theFilamentSpecies[i]->setRadius(VoxelRadius);
    }
  //Normalized lattice voxel radius:
  normLatVoxelRadius = 0.5;
  theVoxelSize = (unsigned int)rint(Length/(VoxelRadius*2));
  startCoord = aStartCoord;
  endCoord = startCoord+Protofilaments*theVoxelSize;
  return endCoord-startCoord;
}

void FilamentProcess::initializeThird()
{
  if(!isCompartmentalized)
    {
      thePoints.resize(Protofilaments*theVoxelSize);
      initializeProtofilaments();
      elongateProtofilaments();
      connectProtofilaments();
      //enlistLatticeVoxels();
      isCompartmentalized = true;
    }
  theVacantSpecies->setIsPopulated();
}

void FilamentProcess::addCompVoxel(unsigned int protoIndex, 
                                   unsigned int voxelIndex, Point& aPoint)
{
  unsigned int aCoord(startCoord+protoIndex*theVoxelSize+voxelIndex);
  Voxel& aVoxel((*theLattice)[aCoord]);
  aVoxel.point = &thePoints[protoIndex*theVoxelSize+voxelIndex];
  *aVoxel.point = aPoint;
  aVoxel.adjoiningCoords = new unsigned int[theAdjoiningCoordSize];
  aVoxel.diffuseSize = 2;
  for(unsigned int i(0); i != theAdjoiningCoordSize; ++i)
    {
      aVoxel.adjoiningCoords[i] = theNullCoord;
    }
  theVacantSpecies->addCompVoxel(aCoord);
}

void FilamentProcess::initializeDirectionVector()
{ 
  /*
   * MEnd = {Mx, My, Mz};(*minus end*) 
   * PEnd = {Px, Py, Pz};(*plus end*)
   * MTAxis = (PEnd - MEnd)/Norm[PEnd - MEnd] (*direction vector along the MT
   * long axis*)
   */
  //Minus end
  M.x = -Length/2;
  M.y = 0;
  M.z = 0;
  //Rotated Minus end
  theSpatiocyteStepper->rotateX(RotateX, &M, -1);
  theSpatiocyteStepper->rotateY(RotateY, &M, -1);
  theSpatiocyteStepper->rotateZ(RotateZ, &M, -1);
  M.x += C.x;
  M.y += C.y;
  M.z += C.z;
  //Direction vector from the Minus end to center
  T.x = C.x-M.x;
  T.y = C.y-M.y;
  T.z = C.z-M.z;
  //Make T a unit vector
  double NormT(sqrt(T.x*T.x+T.y*T.y+T.z*T.z));
  T.x /= NormT;
  T.y /= NormT;
  T.z /= NormT;
  //Rotated Plus end
  P.x = M.x+Length*T.x;
  P.y = M.y+Length*T.y;
  P.z = M.z+Length*T.z;
}

void FilamentProcess::initializeProtofilaments()
{
  initializeDirectionVector();
  Point R; //Initialize a random point on the plane attached at the minus end
  if(M.x != P.x)
    {
      R.y = 10;
      R.z = 30; 
      R.x = (M.x*T.x+M.y*T.y-R.y*T.y+M.z*T.z-R.z*T.z)/T.x;
    }
  else if(M.y != P.y)
    {
      R.x = 10; 
      R.z = 30;
      R.y = (M.x*T.x-R.x*T.x+M.y*T.y+M.z*T.z-R.z*T.z)/T.y;
    }
  else
    {
      R.x = 10; 
      R.y = 30;
      R.z = (M.x*T.x-R.x*T.x+M.y*T.y-R.y*T.y+M.z*T.z)/T.z;
    }
  Point D; //The direction vector from the minus end to the random point, R
  D.x = R.x-M.x;
  D.y = R.y-M.y;
  D.z = R.z-M.z;
  double NormD(sqrt(D.x*D.x+D.y*D.y+D.z*D.z));
  D.x /= NormD;
  D.y /= NormD;
  D.z /= NormD;
  //std::cout << "D.x:" << D.x << " y:" << D.y << " z:" << D.z << std::endl;
  //std::cout << "T.x:" << T.x << " y:" << T.y << " z:" << T.z << std::endl;
  //The start point of the first protofilament:
  Point S; 
  S.x = M.x;
  S.y = M.y;
  S.z = M.z;
  //std::cout << "S.x:" << S.x << " y:" << S.y << " z:" << S.z << std::endl;
  addCompVoxel(0, 0, S);
  for(unsigned int i(1); i != Protofilaments; ++i)
    {
      Point U(S);
      U.x += i*normVoxelRadius*sqrt(3)*D.x;
      U.y += i*normVoxelRadius*sqrt(3)*D.y;
      U.z += i*normVoxelRadius*sqrt(3)*D.z;
      if(i%2 == 1)
        {
          std::cout << "here" << std::endl;
          U.x += normVoxelRadius*T.x;
          U.y += normVoxelRadius*T.y;
          U.z += normVoxelRadius*T.z;
        }
      addCompVoxel(i, 0, U);
    }
}

void FilamentProcess::elongateProtofilaments()
{
  for(unsigned int i(0); i != Protofilaments; ++i)
    {
      Voxel* startVoxel(&(*theLattice)[startCoord+i*theVoxelSize]);
      Point A(*startVoxel->point);
      for(unsigned int j(1); j != theVoxelSize; ++j)
        {
          A.x += (normVoxelRadius*2)*T.x;
          A.y += (normVoxelRadius*2)*T.y;
          A.z += (normVoxelRadius*2)*T.z;
          addCompVoxel(i, j, A);
        }
    }
}

void FilamentProcess::connectProtofilaments()
{
  for(unsigned int i(0); i != theVoxelSize; ++i)
    {
      for(unsigned int j(0); j != Protofilaments; ++j)
        { 
          if(i > 0)
            { 
              connectNorthSouth(i, j);
            }
          else if(Periodic)
            {
              connectPeriodic(j);
            }
          if(j > 0)
            {
              connectEastWest(i, j);
            }
        }
      /*
      if(Protofilaments > 2)
        {
          connectSeamEastWest(i);
        }
        */
      if(i > 0)
        {
          connectNwSw(i);
        }
    }
}

void FilamentProcess::connectPeriodic(unsigned int j)
{
  unsigned int a(startCoord+j*theVoxelSize+theVoxelSize-1);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned int b(startCoord+j*theVoxelSize); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[NORTH] = b;
  adjoin.adjoiningCoords[SOUTH] = a;
  aVoxel.adjoiningSize = 2;
  adjoin.adjoiningSize = 2;
}

void FilamentProcess::connectNorthSouth(unsigned int i, unsigned int j)
{
  unsigned int a(startCoord+j*theVoxelSize+(i-1));
  Voxel& aVoxel((*theLattice)[a]);
  unsigned int b(startCoord+j*theVoxelSize+i);
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[NORTH] = b;
  adjoin.adjoiningCoords[SOUTH] = a;
  aVoxel.adjoiningSize = 2;
  adjoin.adjoiningSize = 2;
}

void FilamentProcess::connectEastWest(unsigned int i, unsigned int j)
{
  unsigned int a(startCoord+j*theVoxelSize+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned int b(startCoord+(j-1)*theVoxelSize+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void FilamentProcess::connectSeamEastWest(unsigned int i)
{
  unsigned int a(startCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned int b(startCoord+(Protofilaments-1)*theVoxelSize+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void FilamentProcess::connectNwSw(unsigned int i)
{
  unsigned int a(startCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned int b(startCoord+(Protofilaments-1)*theVoxelSize+(i-1)); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void FilamentProcess::enlistLatticeVoxels()
{
  for(unsigned int n(startCoord); n != endCoord; ++n)
    {
      Voxel& offVoxel((*theLattice)[n]);
      offVoxel.diffuseSize = offVoxel.adjoiningSize;
      double rA(theSpatiocyteStepper->getMinLatticeSpace());
      if(rA < normVoxelRadius)
        {
         rA = normVoxelRadius;
        } 
      Point center(*offVoxel.point);
      unsigned int aCoord(theSpatiocyteStepper->point2coord(center));
      Point cl(theSpatiocyteStepper->coord2point(aCoord));
      //theSpecies[3]->addMolecule(aVoxel.coord);
      Point bottomLeft(*offVoxel.point);
      Point topRight(*offVoxel.point);
      bottomLeft.x -= rA+center.x-cl.x+theSpatiocyteStepper->getColLength();
      bottomLeft.y -= rA+center.y-cl.y+theSpatiocyteStepper->getLayerLength();
      bottomLeft.z -= rA+center.z-cl.z+theSpatiocyteStepper->getRowLength();
      topRight.x += rA+cl.x-center.x+theSpatiocyteStepper->getColLength()*1.5;
      topRight.y += rA+cl.y-center.y+theSpatiocyteStepper->getLayerLength()*1.5;
      topRight.z += rA+cl.z-center.z+theSpatiocyteStepper->getRowLength()*1.5;
      unsigned int blRow(0);
      unsigned int blLayer(0);
      unsigned int blCol(0);
      theSpatiocyteStepper->point2global(bottomLeft, blRow, blLayer, blCol);
      unsigned int trRow(0);
      unsigned int trLayer(0);
      unsigned int trCol(0);
      theSpatiocyteStepper->point2global(topRight, trRow, trLayer, trCol);
      std::vector<unsigned int> checkedAdjoins;
      for(unsigned int i(blRow); i < trRow; ++i)
        {
          for(unsigned int j(blLayer); j < trLayer; ++j)
            {
              for(unsigned int k(blCol); k < trCol; ++k)
                {
                  unsigned int lat(theSpatiocyteStepper->global2coord(i, j, k));
                  Voxel& latVoxel((*theLattice)[lat]);
                  if(latVoxel.id != theSpatiocyteStepper->getNullID())
                    {
                      //theSpecies[3]->addMolecule(latVoxel);
                      Point aPoint(theSpatiocyteStepper->coord2point(lat));
                      if(inMTCylinder(aPoint))
                        {
                          for(unsigned int l(0); l != theAdjoiningCoordSize;
                              ++l)
                            {
                              unsigned adj(latVoxel.adjoiningCoords[l]);
                              Voxel& adjoin((*theLattice)[adj]);
                              if(adjoin.id == theComp->vacantSpecies->getID())
                                {
                                  checkedAdjoins.push_back(adj);
                                  addDirect(offVoxel, n, adjoin, adj);
                                }
                            }
                        }
                    }
                }
            }
        }
      for(unsigned int i(0); i != checkedAdjoins.size(); ++i)
        {
          (*theLattice)[checkedAdjoins[i]].id = theComp->vacantSpecies->getID();
        }
    }
  for(unsigned int i(0); i != occCoords.size(); ++i)
    {
      Voxel& aVoxel((*theLattice)[occCoords[i]]);
      unsigned int* temp = aVoxel.initAdjoins;
      aVoxel.initAdjoins = aVoxel.adjoiningCoords;
      aVoxel.adjoiningCoords = temp;
      aVoxel.diffuseSize = aVoxel.adjoiningSize;
    }
  for(unsigned int i(0); i != occCoords.size(); ++i)
    {
      Voxel& aVoxel((*theLattice)[occCoords[i]]);
      for(unsigned int i(0); i != aVoxel.adjoiningSize; ++i)
        {
          unsigned int aCoord(aVoxel.adjoiningCoords[i]);
          Voxel& adjoin((*theLattice)[aCoord]);
          if(adjoin.id == theComp->vacantSpecies->getID())
            {
              Point adPoint(theSpatiocyteStepper->coord2point(aCoord));
              if(inMTCylinder(adPoint))
                {
                  std::cout << "error in MT Process" << std::endl;
                }
            }
          else if(adjoin.id != theVacantSpecies->getID() &&
                  adjoin.id != theSpatiocyteStepper->getNullID())
            {
              std::cout << "species error in MT Process" << std::endl;
            }
        }
    }
}

void FilamentProcess::addDirect(Voxel& offVoxel, unsigned a,
                                   Voxel& adjoin, unsigned b)
{
  Point aPoint(*offVoxel.point);
  adjoin.id = tempID;
  Point adPoint(theSpatiocyteStepper->coord2point(b));
  if(!inMTCylinder(adPoint))
    { 
      if(initAdjoins(adjoin)) 
        {
          occCoords.push_back(b);
        }
      double dist(getDistance(&aPoint, &adPoint)); 
      if(dist <= normLatVoxelRadius+normVoxelRadius)
        {
          offVoxel.adjoiningCoords[offVoxel.adjoiningSize++] = b;
          updateAdjoinSize(adjoin);
          adjoin.initAdjoins[adjoin.adjoiningSize++] = a;
        }
      else
        { 
          addIndirect(offVoxel, a, adjoin, b);
        }
    }
}

void FilamentProcess::addIndirect(Voxel& offVoxel, unsigned a,
                                     Voxel& latVoxel, unsigned b)
{
  Point aPoint(*offVoxel.point);
  for(unsigned int i(0); i != theAdjoiningCoordSize; ++i)
    {
      unsigned int aCoord(latVoxel.adjoiningCoords[i]);
      Voxel& adjoin((*theLattice)[aCoord]);
      if(adjoin.id == theComp->vacantSpecies->getID() || adjoin.id == tempID)
        {
          Point adPoint(theSpatiocyteStepper->coord2point(aCoord));
          double dist(getDistance(&aPoint, &adPoint)); 
          if(dist <= normVoxelRadius && inMTCylinder(adPoint))
            { 
              offVoxel.adjoiningCoords[offVoxel.adjoiningSize++] = b; 
              initAdjoins(latVoxel);
              updateAdjoinSize(latVoxel);
              latVoxel.initAdjoins[latVoxel.adjoiningSize++] = a;
            }
        }
    }
}

bool FilamentProcess::initAdjoins(Voxel& aVoxel)
{
  if(aVoxel.initAdjoins == NULL)
    {
      aVoxel.adjoiningSize = 0;
      aVoxel.initAdjoins = new unsigned int[theAdjoiningCoordSize];
      for(unsigned int i(0); i != theAdjoiningCoordSize; ++i)
        {
          unsigned int aCoord(aVoxel.adjoiningCoords[i]);
          Point aPoint(theSpatiocyteStepper->coord2point(aCoord));
          if(!inMTCylinder(aPoint))
            {
              aVoxel.initAdjoins[aVoxel.adjoiningSize++] = aCoord;
            }
        }
      return true;
    }
  return false;
}

void FilamentProcess::updateAdjoinSize(Voxel& aVoxel)
{
 if(aVoxel.adjoiningSize >= theAdjoiningCoordSize)
    {
      unsigned int* temp(new unsigned int[aVoxel.adjoiningSize+1]);
      for(unsigned int i(0); i != aVoxel.adjoiningSize; ++i)
        {
          temp[i] = aVoxel.initAdjoins[i];
        }
      delete[] aVoxel.initAdjoins;
      aVoxel.initAdjoins = temp;
    }
}


bool FilamentProcess::inMTCylinder(Point& N)
{
  /*
  Point E(M);
  Point W(P);
  Point S(M);
  double t((-E.x*N.x-E.y*N.y-E.z*N.z+E.x*S.x+E.y*S.y+E.z*S.z+N.x*W.x-S.x*W.x+N.y*W.y-S.y*W.y+N.z*W.z-S.z*W.z)/(E.x*E.x+E.y*E.y+E.z*E.z-2*E.x*W.x+W.x*W.x-2*E.y*W.y+W.y*W.y-2*E.z*W.z+W.z*W.z));
  double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));
  if(dist < Radius)
    {
      return true;
    }
  */
  return false;
}


// The function returns the result when the point (x,y,z) is rotated about
// the line through (a,b,c) with unit direction vector ⟨u,v,w⟩ by the angle.
void FilamentProcess::rotatePointAlongVector(Point& S, double angle)
{
  double x(S.x);
  double y(S.y);
  double z(S.z);
  double a(M.x);
  double b(M.y);
  double c(M.z);
  double u(T.x);
  double v(T.y);
  double w(T.z);
  double u2(u*u);
  double v2(v*v);
  double w2(w*w);
  double cosT(cos(angle));
  double oneMinusCosT(1-cosT);
  double sinT(sin(angle));
  double xx((a*(v2 + w2) - u*(b*v + c*w - u*x - v*y - w*z)) * oneMinusCosT
                + x*cosT + (-c*v + b*w - w*y + v*z)*sinT);
  double yy((b*(u2 + w2) - v*(a*u + c*w - u*x - v*y - w*z)) * oneMinusCosT
                + y*cosT + (c*u - a*w + w*x - u*z)*sinT);
  double zz((c*(u2 + v2) - w*(a*u + b*v - u*x - v*y - w*z)) * oneMinusCosT
                + z*cosT + (-b*u + a*v - v*x + u*y)*sinT);
  S.x = xx;
  S.y = yy;
  S.z = zz;
}




