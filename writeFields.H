#include <hdf5/hdffile.hpp>
#include <hdf5/hdfdataset.hpp>

#ifndef zcfd_writeFields_H
#define zcfd_writeFields_H
#include "fvMesh.H"


namespace Foam
{

wordHashSet checkzGoodFields
(
    const instant& runTime,
    const fvMesh& mesh,
    const wordList& fieldNames
)
{
    wordHashSet goodFields;
    for (const word& fieldName : fieldNames)
    {
        bool good = false;
        good =
            IOobject
            (
                fieldName,
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false  // no register
            ).typeHeaderOk<volScalarField>(false, false);
        
        reduce(good, andOp<bool>());

        if (good)
        {
            goodFields.insert(fieldName);
        }
    }   
    
    Info<< "good=" << goodFields << nl;
    return goodFields;
}

/*---------------------------------------------------------------------------*\
1. Collect 
\*---------------------------------------------------------------------------*/
word writeVolumeFieldData
(
    const word& velocity_field,
    const word& pressure_field,
    const word& turb_field,
    const instant& runTime,
    const fvMesh& mesh
)
{

    // make a directory called zCFDInterface in the case
    mkDir(mesh.time().rootPath()/mesh.time().caseName()/"zCFDInterface");

    // open a file for the mesh
    std::string filename = (
        mesh.time().rootPath()/
        mesh.time().caseName()/
        "zCFDInterface"/
        mesh.time().caseName() + "_results.h5"
    ).c_str();
    hdf::HDFFile<> file(filename, hdf::HDFFile<>::truncate);

    Info<< "Writing Attributes" << endl;
    boost::shared_ptr<hdf::HDFGroup<> > rung = file.openGroup("run", true);
    std::vector<hsize_t> dims(1, 1);

    rung->createAttribute<int>("Cycle", dims)->writeData(1);
    rung->createAttribute<int>("NumProcessors", dims)->writeData(1);
    rung->createAttribute<int>("TimeCycle", dims)->writeData(0);
    rung->createAttribute<int>("TurbModelSANEG", dims)->writeData(0);
    rung->createAttribute<int>("TurbModelSST", dims)->writeData(0);
    rung->createAttribute<int>("TurbModelSSTTrans", dims)->writeData(0);
    rung->createAttribute<int>("TurbSSTModelLnOmega", dims)->writeData(0);
    rung->createAttribute<int>("OpenFoamExport", dims)->writeData(1);

    Info<<"Writing globalIndex"<<nl;
    std::vector<hsize_t> fileDims(1);
    int numCells = mesh.nCells();
    std::vector<int> cellIndex(numCells,0);
    for(int i = 0; i < numCells;++i){
        cellIndex[i] = i;
    } 
    fileDims[0] = numCells;
    hdf::Slab<1> d1(fileDims);
    rung->createDataset<int>("globalIndex", d1)->writeData(cellIndex);
    
    Info<<"Exporting solution data from Foam"<<nl;
    std::vector<hsize_t> solutionDims(2);
    solutionDims[0] = numCells;
    solutionDims[1] = 5; //number of variables to write
    hdf::Slab<1> d2(solutionDims);
    
    std::vector<double> foam_fields(numCells*5);

    Info<<"Fetching pressure data, variable="<< pressure_field <<nl;
    std::vector<int> pressure(numCells);

    volScalarField fld
    (
        IOobject
        (
                pressure_field,
                runTime.name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
        ),
        mesh
    );
    for (label indx=0; indx<mesh.nCells(); ++indx)
    {
        foam_fields[indx*5] = fld[indx];
        // pressure[indx] = fld[indx];
    }

    Info<<"Fetching turbulance data, variable=" << turb_field <<nl;

    std::vector<double> omega(numCells);

    volScalarField omega_fld
    (
        IOobject
        (
                turb_field,
                runTime.name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
        ),
        mesh
    );
    for (label indx=0; indx<mesh.nCells(); ++indx)
    {
        foam_fields[4+indx*5] = omega_fld[indx];
        // omega[indx] = omega_fld[indx];
    }
    
    Info<<"Fetching velocity data, varible=" << velocity_field  <<nl;

    std::vector<std::vector<double>> velocity(numCells,std::vector<double>(3));

    volVectorField velocity_fld
    (
        IOobject
        (
                velocity_field,
                runTime.name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
        ),
        mesh
    );
    for (label indx=0; indx<mesh.nCells(); ++indx)
    {
        foam_fields[1+indx*5] = velocity_fld[indx].x();
        foam_fields[2+indx*5] = velocity_fld[indx].y();
        foam_fields[3+indx*5] = velocity_fld[indx].z();
        // velocity[indx][0] = velocity_fld[indx].x();
        // velocity[indx][1] = velocity_fld[indx].y();
        // velocity[indx][2] = velocity_fld[indx].z();
    }
    Info<<"Done fetching data"<<nl;

    Info<<"Writing solution data"<<nl;
    rung->createDataset<double>("solution", d2)->writeData(foam_fields);

    Info<<"Done writing solution data"<<nl;
    return filename;
}


} // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif