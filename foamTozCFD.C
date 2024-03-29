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
        "uinf",
        "n",
        "Freestream velocity"
    );
    argList::addOption
    (
        "pinf",
        "n",
        "Freestream Pressure"
    );
    argList::addOption
    (
        "rinf",
        "n",
        "Freestream Density"
    );


    #include "setRootCase.H"
    #include "createTime.H"

    const bool writeResults = args.found("results");
    word velocityField;
    word pressureField;
    word turbField;
    double Vel_inf;
    double P_inf;
    double R_inf;

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
    if (writeResults){
        if (!args.readIfPresent("uinf", Vel_inf))
        {
            Info<<"Error Freestream velocity not set, cannot write results with this. Please supply one with -uinf"<<nl;
            return 1;
        }
        if (!args.readIfPresent("pinf", P_inf))
        {
            Info<<"Freestream pressure not set, assuming sea level and converting freestream velocity. If incorrect supply one with -pinf"<<nl;
            P_inf = 101325.0;
        }
        if (!args.readIfPresent("rinf", R_inf))
        {
            Info<<"Error Freestream density not set, assuming sea level. If incorrect supply one with -rinf"<<nl;
            R_inf = 1.225;
        }
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
        Info<<"Converting results with Ref Velocity="<< Vel_inf << ", Ref pressure=" << P_inf << ", Ref density="<< R_inf<<nl;

        // get the available time-steps
        instantList TimeList = runTime.times();
        label nTimes = TimeList.size();

        for (label n=1; n < nTimes; n++)
        {
            if (TimeList[n].value() > 0)
            {

                Info<< "Time = " << TimeList[n].value() << nl;
                word solution_file = writeVolumeFieldData(velocityField, pressureField, turbField, TimeList[n], mesh, Vel_inf, P_inf, R_inf);
                Info<< "Solution written to "<< solution_file <<nl;
            }
        }
    }
 
    Info<< "End" <<  nl << endl;    
  

    return 0;
}


// ************************************************************************* //
