/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        double timeTemp;   //Declares the temperature at each time as double precision
	const label& nCells = mesh.nCells();                        // Total number of cells in the mesh
        double n1;
	double n2;

	n1 = floor(nCells * 0.33) - 1;      // cell 33 if nCells = 100 (-1 because indexing in C++ starts at 0)
	n2 = floor(nCells * 0.67) - 1;      // cell 67 if nCells = 100 (-1 because indexing in C++ starts at 0)
        timeTemp = runTime.value();     //Grabs the temperature at each time step
        	Info<< "runTimeValue = " << runTime.value() << nl << endl;
	forAll(T,celli)    //OpenFOAM for loop that applies to Temp for all cells starting from 0 to cell_i
        {
           sourceTerm.field()[celli] = 0;
        }
	sourceTerm.field()[n1] = 5e7 * mag(Foam::sin(timeTemp)); //source term in [W*m^3/s]
        sourceTerm.field()[n2] = 5e7 * mag(Foam::sin(timeTemp)); //source term in [W*m^3/s]

	while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                sourceTerm * DT/k   //source term becomes [K/s]
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
