#include <hdf5/hdffile.hpp>
#include <hdf5/hdfdataset.hpp>

#ifndef zcfd_writeFields_H
#define zcfd_writeFields_H
#include "fvMesh.H"


namespace Foam
{

wordHashSet checkGoodFields
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
    
    Info<< "Good=" << goodFields << nl;
    return goodFields;
}

word writeVolumeFieldData
(
    const word& velocity_field,
    const word& pressure_field,
    const word& turb_field,
    const instant& runTime,
    const fvMesh& mesh,
    const double Vel_inf,
    const double P_inf,
    const double R_inf,
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

    int num_vars = 5;
    if(turb_field.size() > 0){
        num_vars = 7;
    }

    std::vector<hsize_t> solutionDims(2);
    solutionDims[0] = numCells;
    solutionDims[1] = num_vars; //number of variables to write
    hdf::Slab<1> d2(solutionDims);
    
    std::vector<double> foam_fields(numCells*num_vars, 1.0);

    Info<<"Fetching pressure data, variable="<< pressure_field <<nl;
    std::vector<int> pressure(numCells);

    // Pressure
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

    // Set density
    for (label indx=0; indx<mesh.nCells(); ++indx)
    {
        foam_fields[0 + indx*num_vars] = 1.0;
    }
    
    Info<<"Fetching velocity data, variable=" << velocity_field  <<nl;

    // Velocity
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
    // Need u_inf and M_inf to non dimensionalise
    double gamma = 1.4;
    auto M_inf = Vel_inf / sqrt(gamma * P_inf/R_inf);

    for (label indx=0; indx<mesh.nCells(); ++indx)
    {
        foam_fields[1+indx*num_vars] = velocity_fld[indx].x()/Vel_inf * sqrt(gamma) * M_inf;
        foam_fields[2+indx*num_vars] = velocity_fld[indx].y()/Vel_inf * sqrt(gamma) * M_inf;
        foam_fields[3+indx*num_vars] = velocity_fld[indx].z()/Vel_inf * sqrt(gamma) * M_inf;

        // velocity[indx][0] = velocity_fld[indx].x();
        // velocity[indx][1] = velocity_fld[indx].y();
        // velocity[indx][2] = velocity_fld[indx].z();
    }

    // Set total energy
    for (label indx=0; indx<mesh.nCells(); ++indx)
    {
        auto p = fld[indx] * R_inf;
        auto u = foam_fields[1+indx*num_vars];
        auto v = foam_fields[2+indx*num_vars];
        auto w = foam_fields[3+indx*num_vars];
        auto rE = p / (gamma - 1.0) + 0.5 * (u*u + v*v + w*w);

        foam_fields[4 + indx*num_vars] = rE;
    }

    if(num_vars > 5){
        Info<<"Fetching turbulence data, variable=k" <<nl;

        volScalarField k_fld
        (
            IOobject
            (
                    "k",
                    runTime.name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
            ),
            mesh
        );
        for (label indx=0; indx<mesh.nCells(); ++indx)
        {
            foam_fields[5+indx*num_vars] = k_fld[indx] / (Vel_inf*Vel_inf) * gamma * M_inf;
        }
        Info<<"Fetching turbulence data, variable=" << turb_field <<nl;

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
            foam_fields[6+indx*num_vars] = omega_fld[indx] * sqrt(P_inf/R_inf);
        }
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