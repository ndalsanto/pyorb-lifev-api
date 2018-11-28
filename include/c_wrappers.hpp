/*
@author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
@date 27-11-2018
*/


namespace PyOrbLifeV
{

extern "C" {
    int solve_parameter( double * _param, FemSpecifics _femSpecifics, bool _computeOnlyDimension = true );

    int build_fom_affine_components( char * _operator, int _q, FemSpecifics _femSpecifics, bool _computeOnlyDimension = true );
}

}
