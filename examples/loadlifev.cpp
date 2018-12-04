#include <iostream>
#include "init_lifev.hpp"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChronoManager.hpp>


int main ( int argc, char** argv )
{
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
    MPI_Comm * my_comm ( new MPI_Comm ( MPI_COMM_WORLD ) );
#else
    Epetra_Comm * my_comm ( new Epetra_SerialComm );
#endif

    PyOrbLifeV::FemSpecifics my_fem_specs;

    my_fem_specs.model = "";
    my_fem_specs.datafile_path = "examples/data";
    my_fem_specs.external_communicator = my_comm;
    my_fem_specs.u = nullptr;
    my_fem_specs.A = nullptr;
    my_fem_specs.f = nullptr;

    PyOrbLifeV::LifeVSimulator my_lifev_simulator;

    my_lifev_simulator.initialize( my_fem_specs );

    // random value of parameter for testing thermal block problem
    double parameter[3] = { 1.5, 2.5, 3.5 };

    my_lifev_simulator.perform_simulation( parameter );

    double residual = my_lifev_simulator.compute_residual( );

    my_lifev_simulator.finalize( );

    return ( 0 );

    if( residual < 1.e-5 )
    {
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return ( 0 );
    }
    else
    {
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return ( EXIT_FAILURE );
    }
}
