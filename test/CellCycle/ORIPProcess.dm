CLASSNAME = 'ORIPProcess'
BASECLASS = 'FluxProcess'
PROPERTIES = []

PROTECTED_AUX = '''
  Connection C0;
  Connection C1;
  Connection P0;
'''

defineMethod( 'initialize', '''
  C0 = getConnection( "C0" );
  C1 = getConnection( "C1" );
  P0 = getConnection( "P0" );
''' )

defineMethod( 'process', '''
  static Real prev( 0 );

  const Real E1( C0.getVariable()->getConcentration() ); // const?
  const Real E2( C1.getVariable()->getConcentration() ); // const?

  Real E( E1 + E2 );

  if(prev > 0.2 && E <= 0.2)
  {
    P0.getVariable()->setValue( 0 );
  }

  prev = E;
''' )
