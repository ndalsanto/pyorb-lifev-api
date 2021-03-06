/*
@author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
@date 27-11-2018
*/

#ifndef __PYORB_LIFEV_C_WRAPPERS__
#define __PYORB_LIFEV_C_WRAPPERS__

#include "init_lifev.hpp"

namespace PyOrbLifeV
{

extern "C" {
    int solve_parameter( double * _param, FemSpecifics _femSpecifics, bool _computeOnlyDimension = true );

    int build_fom_affine_components( char * _operator, int _q, FemSpecifics _femSpecifics, bool _computeOnlyDimension = true );
}

}

#endif
