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

LifeV::Real zero_function( const LifeV::Real& /*t*/, const LifeV::Real& /*x*/, const LifeV::Real& /*y*/,
                const LifeV::Real& /*z*/, const LifeV::ID& /*i*/)
{
    return 0.;
}

LifeV::Real sourceFunction( const LifeV::Real& /*t*/, const LifeV::Real& /*x*/, const LifeV::Real& /*y*/,
                const LifeV::Real& /*z*/, const LifeV::ID& /*i*/)
{
    return 1.;
}

LifeV::Real diffusion( const LifeV::Real& /*t*/, const LifeV::Real& x, const LifeV::Real& y,
                       const LifeV::Real& z, const LifeV::ID& /*i*/)
{
    return 1.;
}

LifeVSimulator::
LifeVSimulator( )
{
}

int
LifeVSimulator::
initialize( FemSpecifics& _femSpecifics )
{
    MPI_Comm * external_mpi_communicator = ( MPI_Comm * ) _femSpecifics.external_communicator;

#ifdef HAVE_MPI
    std::cout << "LifeVSimulator::initialize communicator" << std::endl;
    M_comm.reset( new Epetra_MpiComm ( *( external_mpi_communicator ) ) );
#else
    M_comm.reset( new Epetra_SerialComm );
#endif

    std::cout << "LifeVSimulator::initialize M_dataFile" << std::endl;

    std::string dataFileName( _femSpecifics.datafile_path );
    M_dataFile.reset( new GetPot( dataFileName ) );

    initialize_simulation( );
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
    typedef LifeV::RegionMesh< LifeV::LinearTetra >                 mesh_Type;
    typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra >           FESpace_Type;
    typedef LifeV::ETFESpace< mesh_Type, LifeV::MapEpetra, 3, 1 >   uSpaceETA_Type;

    std::cout << "LifeVSimulator::initialize_simulation reading mesh" << std::endl;

    std::shared_ptr< mesh_Type > fullMeshPtr ( new mesh_Type ( M_comm ) );

    LifeV::MeshData meshData;
    meshData.setup ( *M_dataFile, "mesh");
    readMesh (*fullMeshPtr, meshData);

    LifeV::MeshPartitioner< mesh_Type > meshPart;

    meshPart.doPartition ( fullMeshPtr, M_comm );
    M_localMeshPtr = meshPart.meshPartition();

    // Clearing global mesh
    fullMeshPtr.reset();

    std::cout << "LifeVSimulator::initialize_simulation creating fe spaces" << std::endl;

    // Defining finite elements standard and ET spaces
    M_uFESpace.reset( new FESpace_Type ( M_localMeshPtr, (*M_dataFile)( "finite_element/degree", "P1" ), 1, M_comm ) );
    M_ETuFESpace.reset( new uSpaceETA_Type ( M_localMeshPtr, & ( M_uFESpace->refFE() ), & ( M_uFESpace->fe().geoMap() ), M_comm ) );

}

int
LifeVSimulator::
perform_simulation( double const * const _param )
{
    std::shared_ptr< LifeV::MatrixEpetra< LifeV::Real > > Ah( new LifeV::MatrixEpetra< LifeV::Real >( M_ETuFESpace->map(), 100 ) );
    std::shared_ptr< LifeV::VectorEpetra > fh( new LifeV::VectorEpetra( M_ETuFESpace->map(), LifeV::Repeated ) );





}

}
