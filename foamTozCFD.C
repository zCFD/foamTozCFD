/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Writes out the OpenFOAM mesh in zCFD HDF5 mesh format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "zCFDFvMesh.H"
#include "IOobjectList.H"
#include "volFields.H"

//local files
#include "writeFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write an OpenFOAM mesh in zCFD HDF5 format"
    );
    argList::noParallel();
    argList::addBoolOption
    (
        "results",
        "Export OpenFOAM field data in zCFD HDF5 format"
    );
    argList::addOption
    (
        "v",
        "name",
        "Field name for velocity vector (Default U)"
    );
    argList::addOption
    (
        "p",
        "name",
        "Field name for pressure variable (Default p)"
    );
    argList::addOption
    (
        "turb",
        "name",
        "Field name for turbulance variable (Default omega)"
    );
    argList::addOption
    (
        "vinf",
        "name",
        "Freestream velocity"
    );
    argList::addOption
    (
        "minf",
        "name",
        "Freestream Mach number"
    );


    #include "setRootCase.H"
    #include "createTime.H"

    const bool writeResults = args.found("results");
    word velocityField;
    word pressureField;
    word turbField;

    if (!args.readIfPresent("v", velocityField))
    {
        velocityField = "U";
    }
    if (!args.readIfPresent("p", pressureField))
    {
        pressureField = "p";
    }
    if (!args.readIfPresent("turb", turbField))
    {
        turbField = "omega";
    }

    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    Info<<writeResults<<nl;

    zCFDFvMesh mesh
    (
        IOobject
        (
            zCFDFvMesh::defaultRegion,
            runTime.constant(),
            runTime
        )
    );

    const label& nCells = mesh.nCells();
    const label& nPoints = mesh.nPoints();

    Info<< "Total number of cells = "
        << nCells << nl << endl;

    Info<< "Total number of nodes = "
        << nPoints << nl << endl;

    word mesh_file = mesh.writezCFDMesh();

    Info<< "Mesh written to " << mesh_file;

    if (writeResults)
    {
        Info<< nl << "Writing results in zCFD format"<<nl;
        Info<<"Velocity field name = "<<  velocityField << nl;
        Info<<"Pressure field name = "<<  pressureField << nl;

        // get the available time-steps
        instantList TimeList = runTime.times();
        label nTimes = TimeList.size();

        for (label n=1; n < nTimes; n++)
        {
            if (TimeList[n].value() > 0)
            {

                Info<< "Time = " << TimeList[n].value() << nl;
                word solution_file = writeVolumeFieldData(velocityField, pressureField, turbField, TimeList[n], mesh, Vel_inf, M_inf);
                Info<< "Solution written to "<< solution_file <<nl;
            }
        }
    }
 
    Info<< "End" <<  nl << endl;    
  

    return 0;
}


// ************************************************************************* //
