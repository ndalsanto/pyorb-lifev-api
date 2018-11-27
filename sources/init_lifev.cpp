/*
@author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
@date 27-11-2018
*/

#include "init_lifev.hpp"
#include <lifev/core/LifeV.hpp>
#include <mpi.h>
#include <lifev/core/mesh/MeshData.hpp>

namespace PyOrbLifeV
{

LifeVSimulator::
LifeVSimulator( )
{
}

int
LifeVSimulator::
initialize( FemSpecifics& _femSpecifics )
{
    MPI_Comm * external_mpi_communicator = (MPI_Comm * ) _femSpecifics.external_communicator;
    // char * data_file_path = _femSpecifics.datafile_path;
    // char * _model = _femSpecifics.model;

#ifdef HAVE_MPI
    M_comm.reset( new Epetra_MpiComm ( *( external_mpi_communicator ) ) );
#else
    M_comm.reset( new Epetra_SerialComm );
#endif

    std::string dataFileName( _femSpecifics.datafile_path );
    M_dataFile.reset( new GetPot( dataFileName ) );
}

int
LifeVSimulator::
finalize()
{
}

int
LifeVSimulator::
initialize_simulation( )
{
    std::shared_ptr< LifeV::RegionMesh< LifeV::LinearTetra > > fullMeshPtr ( new LifeV::RegionMesh< LifeV::LinearTetra > ( M_comm ) );

    LifeV::MeshData meshData;
    meshData.setup ( *M_dataFile, "mesh");
    readMesh (*fullMeshPtr, meshData);

    LifeV::MeshPartitioner< LifeV::RegionMesh< LifeV::LinearTetra > > meshPart;

    meshPart.doPartition ( fullMeshPtr, M_comm );
    M_localMeshPtr = meshPart.meshPartition();

    // Clearing global mesh
    fullMeshPtr.reset();
}


}
