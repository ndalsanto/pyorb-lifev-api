/*
@author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
@date 27-11-2018
*/

#include "init_lifev.hpp"
#include <lifev/core/LifeV.hpp>
#include <mpi.h>

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




}
