/*
@author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
@date 27-11-2018
*/
#include "c_wrappers.hpp"

namespace PyOrbLifeV
{

int solve_parameter( double * _param, FemSpecifics _femSpecifics, bool _computeOnlyDimension )
{
    PyOrbLifeV::LifeVSimulator my_lifev_simulator;

    std::cout << "Initializing LifeVSimulator " << std::endl;

    my_lifev_simulator.initialize( _femSpecifics );

    int my_size = my_lifev_simulator.get_fem_dimension( );

    if( _computeOnlyDimension )
        return my_size;

    // random value of parameter for testing thermal block problem
    double parameter[3] = { 1.5, 2.5, 3.5 };

    my_lifev_simulator.perform_simulation( _param );

    double * _u = _femSpecifics.u;

    for( int ii(0); ii < my_size; ++ii )
        _u[ii] = my_lifev_simulator.get_solution( )->epetraVector()[0][ii];

    double norm = my_lifev_simulator.get_solution( )->norm2( );

    std::cout << "Squared norm of solution in solve_parameter is " << norm*norm << std::endl;


    my_lifev_simulator.finalize( );
}

int build_fom_affine_components( char * _operator, int _q, FemSpecifics _femSpecifics, bool _computeOnlyDimension )
{

}

}
