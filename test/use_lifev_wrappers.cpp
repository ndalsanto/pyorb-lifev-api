#include <iostream>
#include "test_function.hpp"
#include "init_lifev.hpp"
#include "c_wrappers.hpp"

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

    my_fem_specs.model = model;
    my_fem_specs.datafile_path = datafile_path;
    my_fem_specs.external_communicator = external_communicator;
    my_fem_specs.u = nullptr;
    my_fem_specs.A = nullptr;
    my_fem_specs.f = nullptr;

    double parameter[3] = { 1.5, 2.5, 3.5 };

    int computed_dimension = PyOrbLifeV::solve_parameter( parameter, my_fem_specs, true );

    std::cout << "Computed dimension is " << computed_dimension << std::endl;

    my_fem_specs.u = new double[computed_dimension];

    PyOrbLifeV::solve_parameter( parameter, my_fem_specs, false );

    double * u = my_fem_specs.u;
    double norm( 0. );

    for( int ii(0); ii < computed_dimension; ++ii )
        norm += u[ii] * u[ii];

    std::cout << "Squared norm of solution is " << norm << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}
