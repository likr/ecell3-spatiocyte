#ifndef __FLUXPROCESS_HPP
#define __FLUXPROCESS_HPP

#include <functional>
#include <algorithm>

#include "libecs/libecs.hpp"
#include "libecs/Process.hpp"
#include "libecs/VariableReference.hpp"
#include "libecs/Variable.hpp"
#include "libecs/Stepper.hpp"

namespace libecs
{

  class FluxProcess         
    :  
    public Process
  {
  
  public:

    FluxProcess()
    {
      ; // do nothing
    }

    ~FluxProcess()
    {
      ; // do nothing
    }

    StringLiteral getClassName() const { return "FluxProcess"; }
    
    virtual void initialize()
    {
      Process::initialize();
    }

    void setFlux( const Real aVelocity )
    {
      setActivity( aVelocity );

      // Increase or decrease variables, skipping zero coefficients.
      std::for_each( theVariableReferenceVector.begin(),
		     theFirstZeroVariableReferenceIterator,
		     std::bind2nd
		     ( std::mem_fun_ref
		       ( &VariableReference::addFlux ), aVelocity ) );

      std::for_each( theFirstPositiveVariableReferenceIterator,
		     theVariableReferenceVector.end(),
		     std::bind2nd
		     ( std::mem_fun_ref
		       ( &VariableReference::addFlux ), aVelocity ) );
    }



  };

}

#endif /* __FluxProcess_HPP */







