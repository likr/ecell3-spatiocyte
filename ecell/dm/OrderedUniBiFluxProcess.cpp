#include "libecs.hpp"
#include "Process.hpp"
#include "Util.hpp"
#include "PropertyInterface.hpp"
#include "PropertySlotMaker.hpp"
#include "System.hpp"
#include "Stepper.hpp"
#include "Variable.hpp"
#include "VariableProxy.hpp"

#include "FluxProcess.hpp"
#include "ecell3_dm.hpp"

#define ECELL3_DM_TYPE Process

USE_LIBECS;

ECELL3_DM_CLASS
  :  
  public FluxProcess
{

  ECELL3_DM_OBJECT;
  
 public:

  ECELL3_DM_CLASSNAME()
    {
      ECELL3_CREATE_PROPERTYSLOT_SET_GET( Real, KcF );
      ECELL3_CREATE_PROPERTYSLOT_SET_GET( Real, KcR );
      ECELL3_CREATE_PROPERTYSLOT_SET_GET( Real, Keq );
      ECELL3_CREATE_PROPERTYSLOT_SET_GET( Real, KmS );
      ECELL3_CREATE_PROPERTYSLOT_SET_GET( Real, KmP0 );
      ECELL3_CREATE_PROPERTYSLOT_SET_GET( Real, KmP1 );
      ECELL3_CREATE_PROPERTYSLOT_SET_GET( Real, KiP );
    }
  
  SIMPLE_SET_GET_METHOD( Real, KcF );
  SIMPLE_SET_GET_METHOD( Real, KcR );
  SIMPLE_SET_GET_METHOD( Real, Keq );
  SIMPLE_SET_GET_METHOD( Real, KmS );
  SIMPLE_SET_GET_METHOD( Real, KmP0 );
  SIMPLE_SET_GET_METHOD( Real, KmP1 );
  SIMPLE_SET_GET_METHOD( Real, KiP );
    
  virtual void initialize()
    {
      FluxProcess::initialize();
      S0 = getVariableReference( "S0" );
      P0 = getVariableReference( "P0" );
      P1 = getVariableReference( "P1" );
      C0 = getVariableReference( "C0" );
    }

  virtual void process()
    {
      Real S0Concentration = S0.getConcentration();
      Real P0Concentration = P0.getConcentration();
      Real P1Concentration = P1.getConcentration();
      
      Real Den( KcR * KmS + KcR * S0Concentration 
		+ KcF * KmP1 * P0Concentration / Keq 
		+ KcF * KmP0 * P1Concentration / Keq 
		+ KcR * S0Concentration * P0Concentration / KiP
		+ KcR * P0Concentration * P1Concentration / Keq );
      Real velocity( KcF * KcR * C0.getValue()
		     * (S0Concentration - P0Concentration * P1Concentration / Keq) / Den );
      setFlux( velocity );
    }

 protected:

  Real KcF;
  Real KcR;
  Real Keq;

  Real KmS;
  Real KmP0;
  Real KmP1;
  Real KiP;

  VariableReference S0;
  VariableReference P0;
  VariableReference P1;
  VariableReference C0;
  
};

ECELL3_DM_INIT;
