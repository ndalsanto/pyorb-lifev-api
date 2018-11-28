#include <iostream>
#include "test_function.hpp"
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
    std::cout << "Entering load lifev example " << std::endl;
    // MPI initialization
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
    MPI_Comm * my_comm ( new MPI_Comm ( MPI_COMM_WORLD ) );
#else
    Epetra_Comm * my_comm ( new Epetra_SerialComm );
#endif

    std::cout << "Declaring FemSpecifics " << std::endl;

    PyOrbLifeV::FemSpecifics my_fem_specs;

    char * model = nullptr;
    char datafile_path[9] = {'t', 'e', 's', 't', '/', 'd', 'a', 't', 'a' };
    MPI_Comm * external_communicator = my_comm;
    double * u = nullptr;
    double * A = nullptr;
    double * f = nullptr;

    my_fem_specs.model = model;
    my_fem_specs.datafile_path = datafile_path;
    my_fem_specs.external_communicator = external_communicator;
    my_fem_specs.u = u;
    my_fem_specs.A = A;
    my_fem_specs.f = f;

    PyOrbLifeV::LifeVSimulator my_lifev_simulator;

    std::cout << "Initializing LifeVSimulator " << std::endl;

    my_lifev_simulator.initialize( my_fem_specs );

    std::cout << "Finalizing LifeVSimulator " << std::endl;

    my_lifev_simulator.finalize( );

    std::cout << "Exiting MPI " << std::endl;

    #ifdef HAVE_MPI
        MPI_Finalize();
    #endif

    return ( EXIT_SUCCESS );
}
