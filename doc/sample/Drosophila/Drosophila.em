
# created by eml2em program
# from file: droTest.eml, date: Sun Oct 13 15:04:54 2002
#

Stepper FixedODE1Stepper( DE )
#Stepper ODE23Stepper( DE )
#Stepper ODE45Stepper( DE )
{
	# no property
}

System CompartmentSystem( / )
{
	StepperID	DE;
	Variable Variable( SIZE ) { Value 0.000000000000001; }
}

System CompartmentSystem( / )
{
	StepperID	DE;

}

System CompartmentSystem( /CELL )
{
	StepperID	DE;

}

System CompartmentSystem( /CELL/CYTOPLASM )
{
	StepperID	DE;

	Variable Variable( SIZE ) { Value 1e-18; }

	Variable Variable( M )
	{
		Value	3.61328202E-01;
	}
	
	Variable Variable( Pn )
	{
		Value	6.21367E-01;
	}
	
	Variable Variable( P0 )
	{
		Value	3.01106835E-01;
	}
	
	Variable Variable( P1 )
	{
		Value	3.01106835E-01;
	}
	
	Variable Variable( P2 )
	{	
		Value	3.61328202E-01;
	}
	
	Process PythonProcess( R_toy1 )
	{
		InitializeMethod "vs = 0.76; KI = 1";
		ProcessMethod "self.setFlux(((vs * KI)/(KI + (C0.Concentration * C0.Concentration * C0.Concentration))) * self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:M 1 ] [ C0 Variable:.:Pn 0 ];
	}
	
	Process PythonProcess( R_toy2 )
	{
		InitializeMethod "vm = 0.65; Km = 0.5";
		ProcessMethod "self.setFlux(((-1 * vm * P0.Concentration)/(Km + P0.Concentration)) * self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:M 1 ];
	}
	
	Process PythonProcess( R_toy3 )
	{
		InitializeMethod "Km = 0.38";
		ProcessMethod "self.setFlux((Km * C0.Concentration) * self.getSuperSystem().SizeN_A)";
		VariableReferenceList	[ P0 Variable:.:P0 1 ] [ C0 Variable:.:M 0 ];
	}
	
	Process PythonProcess( R_toy4 )
	{
		InitializeMethod "V1 = 3.2; K1 = 2";
		ProcessMethod "self.setFlux(((-1 * V1 * C0.Concentration) / (K1 + C0.Concentration)) * self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P0 1 ] [ C0 Variable:.:P0 0 ];
	}
	
	Process PythonProcess( R_toy5 )
	{
		InitializeMethod "V2 = 1.58; K2 = 2";
		ProcessMethod "self.setFlux(((V2 * C0.Concentration) / (K2 + C0.Concentration)) * self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P0 1 ] [ C0 Variable:.:P1 0 ];
	}
	
	Process PythonProcess( R_toy6 )
	{
		InitializeMethod "V1 = 3.2; K1 = 2";
		ProcessMethod "self.setFlux(((V1 * C0.Concentration) / (K1 + C0.Concentration)) * self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P1 1 ] [ C0 Variable:.:P0 0 ];
	}
	
	Process PythonProcess( R_toy7 )
	{
		InitializeMethod "V2 = 1.58; K2 = 2";
		ProcessMethod "self.setFlux(((-1 * V2 * C0.Concentration) / (K2 + C0.Concentration)) * self.getSuperSystem().SizeN_A)";
		VariableReferenceList	[ P0 Variable:.:P1 1 ] [ C0 Variable:.:P1 0 ];
	}
	
	Process PythonProcess( R_toy8 )
	{
		InitializeMethod "V3 = 5; K3 = 2";
		ProcessMethod "self.setFlux(((-1 * V3 * C0.Concentration) / (K3 + C0.Concentration)) * self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P1 1 ] [ C0 Variable:.:P1 0 ];
	}
	
	Process PythonProcess( R_toy9 )
	{
		InitializeMethod "V4 = 2.5; K4 = 2";
		ProcessMethod "self.setFlux(((V4 * C0.Concentration) / (K4 + C0.Concentration)) * self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P1 1 ] [ C0 Variable:.:P2 0 ];
	}
	
	Process PythonProcess( R_toy10 )
	{
		InitializeMethod "V3 = 5; K3 = 2";
		ProcessMethod "self.setFlux(((V3 * C0.Concentration) / (K3 + C0.Concentration)) *  self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P2 1 ] [ C0 Variable:.:P1 0 ];

	}
	
	Process PythonProcess( R_toy11 )
	{
		InitializeMethod "V4 = 2.5; K4 = 2";
		ProcessMethod "self.setFlux(((-1 * V4 * C0.Concentration) / (K4 + C0.Concentration)) *  self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P2 1 ] [ C0 Variable:.:P2 0 ];
	}
	
	Process PythonProcess( R_toy12 )
	{
		InitializeMethod "K1 = 1.9";
		ProcessMethod "self.setFlux((-1 * K1 * C0.Concentration) *  self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P2 1 ] [ C0 Variable:.:P2 0 ];
	}
	
	Process PythonProcess( R_toy13 )
	{
		InitializeMethod "k2 = 1.3";
		ProcessMethod "self.setFlux((k2 * C0.Concentration) *  self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P2 1 ] [ C0 Variable:.:Pn 0 ];
	}
	
	Process PythonProcess( R_toy14 )
	{
		InitializeMethod "vd = 0.95; Kd = 0.2";
		ProcessMethod "self.setFlux(((-1 * vd * C0.Concentration) / (Kd + C0.Concentration)) *  self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:P2 1 ] [ C0 Variable:.:P2 0 ];
	}
	
	Process PythonProcess( R_toy15 )
	{
		InitializeMethod "k1 = 1.9";
		ProcessMethod "self.setFlux((k1 * C0.Concentration) *  self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:Pn 1 ] [ C0 Variable:.:P2 0 ];
	}
	
	Process PythonProcess( R_toy16 )
	{
		InitializeMethod "k2 = 1.3";
		ProcessMethod "self.setFlux((-1 * k2 * C0.Concentration) *  self.getSuperSystem().SizeN_A)";

		VariableReferenceList	[ P0 Variable:.:Pn 1 ] [ C0 Variable:.:Pn 0 ];
	}
	
	
}

