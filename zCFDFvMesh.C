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

\*---------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <vector>
#include <hdf5/hdffile.hpp>
#include <hdf5/hdfdataset.hpp>

using std::ofstream;
using std::ios;

#include "Time.H"
#include "zCFDFvMesh.H"
#include "primitiveMesh.H"
#include "wallFvPatch.H"
#include "symmetryFvPatch.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zCFDFvMesh::zCFDFvMesh(const IOobject& io)
:
    fvMesh(io)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zCFDFvMesh::writezCFDMesh() const
{
    // make a directory called proInterface in the case
    mkDir(time().rootPath()/time().caseName()/"zCFDInterface");

    // open a file for the mesh
    std::string filename = (
        time().rootPath()/
        time().caseName()/
        "zCFDInterface"/
        time().caseName() + ".h5"
    ).c_str();
    hdf::HDFFile<> file(filename, hdf::HDFFile<>::truncate);

    Info<< "Writing Header" << endl;

    boost::shared_ptr<hdf::HDFGroup<> > meshg = file.openGroup("mesh", true);

    std::vector<hsize_t> dims(1, 1);
    meshg->createAttribute<int>("numFaces", dims)->writeData(nFaces());
    meshg->createAttribute<int>("numCells", dims)->writeData(nCells());

    int numFaces = nFaces();
    int numCells = nCells();

    std::vector<hsize_t> fileDims(2);

   Info<< "Writing Points" << endl;


    // Writing points
    {
      std::vector<double> coord;
      coord.reserve(nPoints()*3);

      const pointField& p = points();

       forAll(p, pointI)
       {
         coord.push_back(p[pointI].x());
         coord.push_back(p[pointI].y());
         coord.push_back(p[pointI].z());
       }

      fileDims[0] = coord.size()/3;
      fileDims[1] = 3;
      hdf::Slab<1> d1(fileDims);
      meshg->createDataset<double>("nodeVertex", d1)->writeData(coord);
    }

       Info<< "Writing Neighbours" << endl;


    const labelUList& own = owner();
    const labelUList& nei = neighbour();


    //const labelUList& cellZone = cellZones();
/*
    const wordList& cellZone = cellZones().names();
    forAll(cellZone, cz)
    {
      Info << cellZone[cz] << endl;

      Info << cellZones().whichZone(0) << endl;
      //Info<< cellZones()[cz] << std::endl;
    }
*/
    std::vector<int> cellZone(numCells,0);
    
    //std::map<int,int> foamCellZoneMap;
    for(int i = 0; i < numCells;++i){
       int cz = cellZones().whichZone(i);

        cellZone[i] = cz+1;
    } 
{
    fileDims[0] = cellZone.size();
    fileDims[1] = 1;
    hdf::Slab<1> d1(fileDims);
    meshg->createDataset<int>("cellZone", d1)->writeData(cellZone);
}

    Info<< "Face Zones" << endl;
    const wordList& faceZone = faceZones().names();
    forAll(faceZone, cz)
    {
      Info << faceZone[cz] << endl;
      //Info<< cellZones()[cz] << std::endl;
      //    }
    }

    { // Face cell
      fileDims[0] = numFaces;
      fileDims[1] = 2;

      hdf::Slab<1> d2(fileDims);

      dims[0] = 2 * numFaces;
      hdf::Slab<1> memFaceCell(dims);

      boost::shared_ptr<hdf::HDFDataSet<> > dataset = meshg->createDataset<int>("faceCell", d2);

      std::vector<int> faceCell;
      //faceCell.reserve(2 * 50000);

      int faceOffset=0;
      int i=0;
      forAll(own, faceI)
      {
        int left = own[faceI];
        int right = nei[faceI];

        faceCell.push_back(left);
        faceCell.push_back(right);
/*
        if(i % 50000)
        {
          std::vector<hsize_t> dims(2,0);
          dims[0] = numFaces;
          dims[1] = 2;
          hdf::Slab<2> filespace(dims);
          dims[0] = 50000;
          dims[1] = 2;
          hdf::Slab<2> memspace(dims);
          std::vector<hsize_t> offset(2,0);
          offset[0] = faceOffset;
          std::vector<hsize_t> stride(2,1);
          hdf::Slab<2> fileselec(filespace,offset,stride,dims);
          dataset->writeData(&faceCell[0],memspace,fileselec);

          faceCell.resize(0);
          faceOffset += 50000;
        }*/
       i++;
      }

      forAll(boundary(), patchI)
      {
          const faceUList& patchFaces = boundaryMesh()[patchI];

          const labelList& patchFaceCells =
              boundaryMesh()[patchI].faceCells();
          forAll(patchFaces, faceI)
          {
            int left = patchFaceCells[faceI];
            int right = numCells;

            faceCell.push_back(left);
            faceCell.push_back(right);

            numCells++;
          }
      }
      if(faceCell.size())
      {
        std::vector<hsize_t> dims(2,0);
        dims[0] = numFaces;
        dims[1] = 2;
        hdf::Slab<2> filespace(dims);
        dims[0] = faceCell.size()/2;
        dims[1] = 2;
        hdf::Slab<2> memspace(dims);
        std::vector<hsize_t> offset(2,0);
        offset[0] = faceOffset;
        std::vector<hsize_t> stride(2,1);
        hdf::Slab<2> fileselec(filespace,offset,stride,dims);
        dataset->writeData(&faceCell[0],memspace,fileselec);
      }
      faceOffset += faceCell.size()/2;
      faceCell.resize(0);

    }

    Info<< "Writing Faces" << endl;

    const faceList& fcs = faces();

    size_t totalFaceNodes=0;
    { // FaceType
      std::vector<int> faceType;
      faceType.reserve(numFaces);
      forAll(own, faceI)
      {
        const labelList& l = fcs[faceI];
        int n = l.size();
        faceType.push_back(n);
        totalFaceNodes+=n;
      }
      forAll(boundary(), patchI)
      {
          const faceUList& patchFaces = boundaryMesh()[patchI];
          forAll(patchFaces, faceI)
          {
              const labelList& l = patchFaces[faceI];
              int n = l.size();
              faceType.push_back(n);
              totalFaceNodes+=n;
          }
      }
      fileDims[0] = numFaces;
      fileDims[1] = 1;
      hdf::Slab<1> d1(fileDims);
      meshg->createDataset<int>("faceType", d1)->writeData(faceType);
    }
    { // Face nodes
      fileDims[0] = totalFaceNodes;
      fileDims[1] = 1;
      hdf::Slab<1> d1(fileDims);
      boost::shared_ptr<hdf::HDFDataSet<> > dataset = meshg->createDataset<int>("faceNodes", d1);

      std::vector<int> faceNodes;
      faceNodes.reserve(totalFaceNodes);
      forAll(own, faceI)
      {
        const labelList& l = fcs[faceI];
        forAll(l, lI)
        {
          faceNodes.push_back(l[lI]);
        }
      }
      forAll(boundary(), patchI)
      {
          const faceUList& patchFaces = boundaryMesh()[patchI];
          forAll(patchFaces, faceI)
          {
              const labelList& l = patchFaces[faceI];
              forAll(l, lI)
              {
                faceNodes.push_back(l[lI]);
              }
          }
      }
      dataset->writeData(faceNodes);
    }

    Info<< "Writing Boundary Conditions" << endl;

    { // FaceBC
      std::vector<int> faceBC(numFaces, 0);
      int i=0;
      forAll(own, faceI)
      {
        faceBC[i] = 0; i++;
      }
      forAll(boundary(), patchI)
      {
          const faceUList& patchFaces = boundaryMesh()[patchI];
          forAll(patchFaces, faceI)
          {
            int bc = -1;
            // Write patch type
             if (isA<wallFvPatch>(boundary()[patchI]))
             {
                 bc=3;
             }
             else if (isA<symmetryFvPatch>(boundary()[patchI]))
             {
                 bc=7;
             }
             else
             {
                 bc=9;
             }

            faceBC[i] = bc; i++;
          }
      }
      assert(i == numFaces);
      fileDims[0] = numFaces;
      fileDims[1] = 1;
      hdf::Slab<1> d1(fileDims);
      meshg->createDataset<int>("faceBC", d1)->writeData(faceBC);
    }
    { // Face Info
      std::vector<int> faceInfo;
      faceInfo.reserve(2 * numFaces);
      forAll(own, faceI)
      {
        int z = 0;
        int dummy = 0;
        faceInfo.push_back(z);
        faceInfo.push_back(dummy);
      }
      forAll(boundary(), patchI)
      {
          const faceUList& patchFaces = boundaryMesh()[patchI];
          forAll(patchFaces, faceI)
          {

        int z = patchI+1;
        int dummy = 0;
        faceInfo.push_back(z);
        faceInfo.push_back(dummy);
         }
      }
      fileDims[0] = numFaces;
      fileDims[1] = 2;
      hdf::Slab<1> d1(fileDims);
      meshg->createDataset<int>("faceInfo", d1)->writeData(faceInfo);
    }

}


// ************************************************************************* //
