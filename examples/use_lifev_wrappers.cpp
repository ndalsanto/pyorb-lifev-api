#include <iostream>
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
    // MPI initialization
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

    double parameter[3] = { 1.5, 2.5, 3.5 };

    int computed_dimension = PyOrbLifeV::solve_parameter( parameter, my_fem_specs, true );

    my_fem_specs.u = new double[computed_dimension];

    PyOrbLifeV::solve_parameter( parameter, my_fem_specs, false );

    double * u = my_fem_specs.u;
    double norm2( 0. );

    for( int ii(0); ii < computed_dimension; ++ii )
        norm2 += u[ii] * u[ii];

    std::cout << "Squared norm of solution is " << norm2 << std::endl;

    // check on the norm of the solution
    if( std::abs( norm2 - 0.784939 ) < 1.e-4 )
    {
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return ( EXIT_SUCCESS );
    }
    else
    {
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return ( EXIT_FAILURE );
    }
}
