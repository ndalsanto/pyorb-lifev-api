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
    PyOrbLifeV::LifeVSimulator my_lifev_simulator;

    std::cout << "Initializing LifeVSimulator " << std::endl;

    my_lifev_simulator.initialize( _femSpecifics );

    if( *_operator == 'A' )
    {
        double * affine_component_q = _femSpecifics.A;
        double * fake_parameter = new double[3];

        fake_parameter[0] = 0.0; fake_parameter[1] = 0.0; fake_parameter[2] = 0.0;

        if( _q < 3 )
            fake_parameter[_q] = 1.0;

        my_lifev_simulator.build_stiffness_matrix( fake_parameter );

        delete [] fake_parameter;

        auto AqPointer = my_lifev_simulator.get_stiffness( )->matrixPtr( );

        if( _computeOnlyDimension )
            return AqPointer->NumMyNonzeros( );

        double * values;
        int * indices;
        int counter = 0;

        for( int iR(0); iR < AqPointer->NumMyRows( ); ++iR )
        {
            int row_entries = AqPointer->NumMyEntries(iR);
            int num_entries = 0;

            // the number of nonzeros per row changes with every row.
            values  = new double[row_entries];
            indices = new int[row_entries];

            AqPointer->ExtractMyRowCopy( iR, row_entries, num_entries, values, indices );

            for( int iV(0); iV < row_entries; ++iV )
            {
                affine_component_q[counter] = AqPointer->GRID( iR );
                counter++;
                affine_component_q[counter] = AqPointer->GCID( indices[iV] );
                counter++;
                affine_component_q[counter] = values[iV];
                counter++;
            }

            delete [] values;
            delete [] indices;
        }
    }
    else if( *_operator == 'f' )
    {
        double * affine_component_fq = _femSpecifics.f;

        int my_size = my_lifev_simulator.get_fem_dimension( );

        if( _computeOnlyDimension )
            return my_size;

        double * fake_parameter = new double[3];
        fake_parameter[0] = 0.0; fake_parameter[1] = 0.0; fake_parameter[2] = 0.0;
        my_lifev_simulator.build_fem_vector( fake_parameter );

        auto fqPointer = my_lifev_simulator.get_rhs( );

        for( int ii(0); ii < my_size; ++ii )
            affine_component_fq[ii] = fqPointer->epetraVector()[0][ii];
    }

    my_lifev_simulator.finalize( );
}

}
