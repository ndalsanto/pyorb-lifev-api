#include <iostream>
#include "test_function.hpp"
#include "init_lifev.hpp"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChronoManager.hpp>


int main (int argc, char **argv)
{
    std::cout << "Hello world!" << std::endl;

    LifeV::LifeChrono globalChrono;

    test_function( 3 );

    PyOrbLifeV::FemSpecifics my_fem_specs;

    char * model;
    char * datafile_path;
    void * external_communicator;
    double * u;
    double * A;
    double * f;

    my_fem_specs.model = model;
    my_fem_specs.datafile_path = datafile_path;
    my_fem_specs.external_communicator = external_communicator;
    my_fem_specs.u = u;
    my_fem_specs.A = A;
    my_fem_specs.f = f;

    PyOrbLifeV::LifeVSimulator my_lifev_simulator;

    my_lifev_simulator.initialize( my_fem_specs );

    my_lifev_simulator.finalize( );

    return 0;
}
