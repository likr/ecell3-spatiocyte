#!/usr/bin/python2

from ecssupport import *

aSimulator = aMainWindow.theSimulator

aMainWindow.printMessage("Load Rule...\n")

aSimulator.createEntity( 'System', ( SYSTEM, '/', 'CELL' ), 'cell' )
aSimulator.createEntity( 'System', ( SYSTEM, '/CELL', 'CYTOPLASM' ), 'cytoplasm' )
aSimulator.createEntity( 'System', ( SYSTEM, '/CELL', 'MEMBRANE' ), 'membrane' )
aSimulator.createEntity( 'System', ( SYSTEM, '/', 'ENVIRONMENT' ), 'environment' )

aMainWindow.printMessage('make Substances...')

aSimulator.createEntity( 'Substance', ( SUBSTANCE, '/CELL/CYTOPLASM', 'ATP' ), 'substance ATP' )
aSimulator.createEntity( 'Substance', ( SUBSTANCE, '/CELL/CYTOPLASM', 'ADP' ), 'substance ADP' )
aSimulator.createEntity( 'Substance', ( SUBSTANCE, '/CELL/CYTOPLASM', 'AMP' ), 'substance AMP' )

aSimulator.createEntity( 'Substance', ( SUBSTANCE, '/ENVIRONMENT', 'GLU' ), 'Glucose' )
aSimulator.createEntity( 'Substance', ( SUBSTANCE, '/ENVIRONMENT', 'PYR' ), 'Pyruvate' )
aSimulator.createEntity( 'Substance', ( SUBSTANCE, '/ENVIRONMENT', 'LCT' ), 'Lactate' )

aSimulator.createEntity( 'Substance', ( SUBSTANCE, '/CELL/MEMBRANE', 'CI' ), 'Channel 1' )
aSimulator.createEntity( 'Substance', ( SUBSTANCE, '/CELL/MEMBRANE', 'CII' ), 'Channel 2' )
aSimulator.createEntity( 'Substance', ( SUBSTANCE, '/CELL/MEMBRANE', 'CIII' ), 'Channel 3' )


aMainWindow.printMessage("make Reactors...\n")
try:
    aSimulator.createEntity('ConstantActivityReactor',
                                         ( REACTOR, '/CELL/CYTOPLASM', 'RC1' ),
                                         'constant reactor' )
except:
   aMainWindow.printMessage("cannot instantiate ConstantActivityReactor\n")
 
aMainWindow.printMessage("set Quantity...\n")

aSimulator.setProperty( ( SUBSTANCE, '/CELL/CYTOPLASM', 'ATP', 'Quantity' ), (100,) )
aSimulator.setProperty( ( SUBSTANCE, '/CELL/CYTOPLASM', 'ADP' , 'Quantity' ), (120,) )
aSimulator.setProperty( ( SUBSTANCE, '/CELL/CYTOPLASM', 'AMP', 'Quantity' ), (140,) )

aSimulator.setProperty( ( SUBSTANCE, '/ENVIRONMENT', 'GLU', 'Quantity' ), (300,) )
aSimulator.setProperty( ( SUBSTANCE, '/ENVIRONMENT', 'PYR', 'Quantity' ), (400,) )
aSimulator.setProperty( ( SUBSTANCE, '/ENVIRONMENT', 'LCT', 'Quantity' ), (500,) )

aSimulator.setProperty( ( SUBSTANCE, '/CELL/MEMBRANE', 'CI', 'Quantity' ), (10,) )
aSimulator.setProperty( ( SUBSTANCE, '/CELL/MEMBRANE', 'CII', 'Quantity' ), (20,) )
aSimulator.setProperty( ( SUBSTANCE, '/CELL/MEMBRANE', 'CIII', 'Quantity' ), (30,) )

aMainWindow.printMessage("initialize()...\n")
aSimulator.initialize()

#  aMainWindow.printAllProperties( ( SYSTEM, '', '/' ) )
#  aMainWindow.printAllProperties( ( SYSTEM, '/', 'CYTOPLASM' ) )

#  aMainWindow.printProperty( ( SUBSTANCE, '/', 'A', 'Quantity' ) )

aMainWindow.printProperty( ( SUBSTANCE, '/CELL/CYTOPLASM', 'ATP', 'Quantity' ) )
aMainWindow.printProperty( ( SUBSTANCE, '/CELL/CYTOPLASM', 'ADP', 'Quantity' ) )
aMainWindow.printProperty( ( SUBSTANCE, '/CELL/CYTOPLASM', 'AMP', 'Quantity' ) )

aMainWindow.printProperty( ( SUBSTANCE, '/ENVIRONMENT', 'GLU', 'Quantity' ) )
aMainWindow.printProperty( ( SUBSTANCE, '/ENVIRONMENT', 'PYR', 'Quantity' ) )
aMainWindow.printProperty( ( SUBSTANCE, '/ENVIRONMENT', 'LCT', 'Quantity' ) )

aMainWindow.printProperty( ( SUBSTANCE, '/CELL/MEMBRANE', 'CI', 'Quantity' ) )
aMainWindow.printProperty( ( SUBSTANCE, '/CELL/MEMBRANE', 'CII', 'Quantity' ) )
aMainWindow.printProperty( ( SUBSTANCE, '/CELL/MEMBRANE', 'CIII', 'Quantity' ) )







