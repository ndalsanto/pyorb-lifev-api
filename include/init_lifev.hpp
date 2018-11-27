/*
@author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
@date 27-11-2018
*/

#ifndef __INIT_LIFEV_PYORB_API__
#define __INIT_LIFEV_PYORB_API__

#include <iostream>
#include <memory>
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

namespace PyOrbLifeV
{

struct FemSpecifics
{
//public:
      char * model;
      char * datafile_path;
      void * external_communicator;
      double * u;
      double * A;
      double * f;
};

class LifeVSimulator
{

public:

    LifeVSimulator( );

    ~LifeVSimulator( ) { };

    // methods to initialize/finalize LifeV
    int initialize( FemSpecifics& _femSpecifics );

    int finalize();

    // methods to run a LifeV simulation
    int initialize_simulation( );

    // members to initialize LifeV
    std::shared_ptr< Epetra_Comm >                      M_comm;

    // members to initialize a LifeV simulation
    std::shared_ptr<GetPot>                                           M_dataFile;
    std::shared_ptr< LifeV::RegionMesh< LifeV::LinearTetra > >        M_localMeshPtr;
};

}



#endif
