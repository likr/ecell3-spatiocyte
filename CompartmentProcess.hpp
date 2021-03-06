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


#ifndef __CompartmentProcess_hpp
#define __CompartmentProcess_hpp

#include <sstream>
#include <MethodProxy.hpp>
#include "SpatiocyteProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_CLASS(CompartmentProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(CompartmentProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Integer, Filaments);
      PROPERTYSLOT_SET_GET(Integer, Periodic);
      PROPERTYSLOT_SET_GET(Integer, Subunits);
      PROPERTYSLOT_SET_GET(Real, DiffuseRadius);
      PROPERTYSLOT_SET_GET(Real, Length);
      PROPERTYSLOT_SET_GET(Real, LipidRadius);
      PROPERTYSLOT_SET_GET(Real, OriginX);
      PROPERTYSLOT_SET_GET(Real, OriginY);
      PROPERTYSLOT_SET_GET(Real, OriginZ);
      PROPERTYSLOT_SET_GET(Real, RotateX);
      PROPERTYSLOT_SET_GET(Real, RotateY);
      PROPERTYSLOT_SET_GET(Real, RotateZ);
      PROPERTYSLOT_SET_GET(Real, SubunitRadius);
      PROPERTYSLOT_SET_GET(Real, Width);
    }
  CompartmentProcess():
    isCompartmentalized(false),
    dimension(1),
    Filaments(1),
    Periodic(0),
    Subunits(1),
    DiffuseRadius(0),
    Length(0),
    LipidRadius(0),
    nVoxelRadius(0.5),
    OriginX(0),
    OriginY(0),
    OriginZ(0),
    RotateX(0),
    RotateY(0),
    RotateZ(0),
    SubunitRadius(0),
    Width(0),
    theLipidSpecies(NULL),
    theVacantSpecies(NULL) {}
  virtual ~CompartmentProcess() {}
  SIMPLE_SET_GET_METHOD(Integer, Filaments);
  SIMPLE_SET_GET_METHOD(Integer, Periodic);
  SIMPLE_SET_GET_METHOD(Integer, Subunits);
  SIMPLE_SET_GET_METHOD(Real, DiffuseRadius);
  SIMPLE_SET_GET_METHOD(Real, Length);
  SIMPLE_SET_GET_METHOD(Real, LipidRadius);
  SIMPLE_SET_GET_METHOD(Real, OriginX);
  SIMPLE_SET_GET_METHOD(Real, OriginY);
  SIMPLE_SET_GET_METHOD(Real, OriginZ);
  SIMPLE_SET_GET_METHOD(Real, RotateX);
  SIMPLE_SET_GET_METHOD(Real, RotateY);
  SIMPLE_SET_GET_METHOD(Real, RotateZ);
  SIMPLE_SET_GET_METHOD(Real, SubunitRadius);
  SIMPLE_SET_GET_METHOD(Real, Width);
  virtual void prepreinitialize()
    {
      SpatiocyteProcess::prepreinitialize();
    }
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      theInterfaceVariable = createVariable("Interface");
      theVacantVariable = createVariable("Vacant");
      if(LipidRadius)
        {
          if(LipidRadius < 0)
            {
              LipidRadius = 0;
            }
          else
            {
              theLipidVariable = createVariable("Lipid");
            }
        }
      SpatiocyteProcess::initialize();
      theInterfaceSpecies = theSpatiocyteStepper->addSpecies(
                                                       theInterfaceVariable);
      theVacantSpecies = theSpatiocyteStepper->addSpecies(theVacantVariable);
      if(LipidRadius)
        {
          theLipidSpecies = theSpatiocyteStepper->addSpecies(
                                                       theLipidVariable);
        }
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          Species* aSpecies(theSpatiocyteStepper->variable2species(
                                   (*i).getVariable())); 
          if((*i).getCoefficient())
            {
              theLipidCompSpecies.push_back(aSpecies);
            }
          else
            {
              theVacantCompSpecies.push_back(aSpecies);
            }
        }
      if(!DiffuseRadius)
        {
          DiffuseRadius = theSpatiocyteStepper->getVoxelRadius();
        }
      if(!SubunitRadius)
        {
          SubunitRadius = DiffuseRadius;
        }
      //Lattice voxel radius:
      VoxelRadius = theSpatiocyteStepper->getVoxelRadius();
      //Normalized off-lattice voxel radius:
      nSubunitRadius = SubunitRadius/(VoxelRadius*2);
      nDiffuseRadius = DiffuseRadius/(VoxelRadius*2);
      //Normalized lipid voxel radius:
      nLipidRadius = LipidRadius/(VoxelRadius*2);
    }
  virtual void initializeFirst()
    {
      theComp = new Comp;
      theVacantSpecies->setIsCompVacant();
      theVacantSpecies->setIsOffLattice();
      theVacantSpecies->setComp(theComp);
      if(theLipidSpecies)
        {
          theLipidSpecies->setIsCompVacant();
          theLipidSpecies->setIsOffLattice();
          theLipidSpecies->setComp(theComp);
        }
      for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
        {
          theLipidCompSpecies[i]->setIsOffLattice();
          //setVacantSpecies must be declared here since it needs
          //to be overwritten by DiffusionProcess in initializeSecond:
          theLipidCompSpecies[i]->setVacantSpecies(theLipidSpecies);
          theLipidCompSpecies[i]->setComp(theComp);
        }
      for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
        {
          theVacantCompSpecies[i]->setIsOffLattice();
          //setVacantSpecies must be declared here since it needs
          //to be overwritten by DiffusionProcess in initializeSecond:
          theVacantCompSpecies[i]->setVacantSpecies(theVacantSpecies);
          theVacantCompSpecies[i]->setComp(theComp);
          if(theLipidSpecies)
            {
              theVacantCompSpecies[i]->setIsMultiscale();
            }
        }
    }
  virtual unsigned getLatticeResizeMol(unsigned);
  virtual void initializeThird();
  virtual void printParameters();
  virtual void updateResizedLattice();
  void addCompVoxel(unsigned, unsigned, Point&, Species*, unsigned, unsigned);
  void initializeVectors();
  void initializeFilaments(Point&, unsigned, unsigned, double, Species*,
                           unsigned);
  void elongateFilaments(Species*, unsigned, unsigned, unsigned, double);
  void connectSubunit(unsigned, unsigned);
  void connectFilaments(unsigned, unsigned, unsigned);
  void addInterfaceVoxel(unsigned, unsigned);
  void setCompartmentDimension();
  void setVacantCompSpeciesProperties();
  void setVacantCompMultiscaleProperties();
  void setLipidCompSpeciesProperties();
  void setDiffuseSize(unsigned, unsigned);
  void interfaceSubunits();
  void enlistInterfaceVoxels();
  void enlistNonIntersectInterfaceVoxels();
  void enlistSubunitInterfaceAdjoins();
  void addNonIntersectInterfaceVoxel(unsigned, Point&);
  void rotate(Point&);
  void addAdjoin(unsigned, unsigned);
  bool isInside(Point&);
  void setSpeciesIntersectLipids();
  Point getStartVoxelPoint();
  int getCoefficient(Species*);
  Species* coefficient2species(int);
protected:
  bool isCompartmentalized;
  unsigned dimension;
  unsigned endMol;
  unsigned Filaments;
  unsigned Periodic;
  unsigned LipidRows;
  unsigned LipidCols;
  unsigned lipStartMol;
  unsigned subStartMol;
  unsigned Subunits;
  double DiffuseRadius;
  double Height;
  double nHeight;
  double Length;
  double lengthDisplace;
  double lengthDisplaceOpp;
  double LipidRadius;
  double nDiffuseRadius;
  double nLength;
  double nLipidRadius;
  double nSubunitRadius;
  double nVoxelRadius;
  double nWidth;
  double OriginX;
  double OriginY;
  double OriginZ;
  double RotateX;
  double RotateY;
  double RotateZ;
  double SubunitRadius;
  double surfaceDisplace;
  double VoxelRadius;
  double Width;
  double widthDisplace;
  double widthDisplaceOpp;
  Comp* theComp;
  Point heightEnd;
  Point heightVector;
  Point lengthEnd;
  Point lengthStart;
  Point lengthVector;
  Point lipidStart;
  Point Origin;
  Point subunitStart;
  Point surfaceNormal;
  Point widthEnd;
  Point widthVector;
  Species* theLipidSpecies;
  Species* theInterfaceSpecies;
  Species* theVacantSpecies;
  Variable* theInterfaceVariable;
  Variable* theLipidVariable;
  Variable* theVacantVariable;
  std::vector<Point> thePoints;
  std::vector<Species*> theLipidCompSpecies;
  std::vector<Species*> theVacantCompSpecies;
  std::vector<std::vector<unsigned> > subunitInterfaces;
};

#endif /* __CompartmentProcess_hpp */




