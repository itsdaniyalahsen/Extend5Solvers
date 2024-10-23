/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    acousticFoam

Description
    Fully Implict Solver for the Acoustic Pressure Equation. 
    Uses SIMPLE control to solve the acoustic pressure equation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    pimpleControl pimple(mesh);

#   include "createFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating solution using Wave Equation Solver\n" << endl;

//#   include "CourantNo.H"
//
        while (runTime.run())
        {

//#           include "readTimeControls.H"
//#           include "setDeltaT.H"

            runTime++;

            Info << "Time = " << runTime.timeName() << nl << endl;

            // pimple loop
            while (pimple.loop()) {
              // piso loop
                while(pimple.correct())
                {
                    fvScalarMatrix pEqn
                    (
                      fvm::d2dt2(p) - sqr(c) * fvm::laplacian(p) + nu * mag(fvc::grad(p))
                    );

                    pEqn.relax();
                    pEqn.solve();
                }
            }

        runTime.write();

        Info << "Execution Time = " << runTime.elapsedCpuTime() << " s"
          << " Clock Time = " << runTime.elapsedClockTime() << " s" 
          << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
