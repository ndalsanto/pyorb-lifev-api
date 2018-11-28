#ifndef _CUSTOM_FUNCTOR_HPP_
#define _CUSTOM_FUNCTOR_HPP_

#include <lifev/core/LifeV.hpp>

namespace PyOrbLifeV
{

template < typename Return_Type>
class customFunctor
{

public:
    typedef Return_Type return_Type;

    return_Type operator() ( const LifeV::VectorSmall<3> spaceCoordinates )
    {
        LifeV::Real x = spaceCoordinates[0];
        LifeV::Real y = spaceCoordinates[1];
        LifeV::Real z = spaceCoordinates[2];

        return myFunction( 0, x, y, z, 0 );
    }

    customFunctor ( return_Type (*newFunction)(const LifeV::Real& /*t*/, const LifeV::Real& x, const LifeV::Real& y,
                                               const LifeV::Real& z, const LifeV::ID& /*i*/) )
        : myFunction( newFunction ) {  }
    ~customFunctor() {}

private:

    return_Type (*myFunction)( const LifeV::Real& /*t*/, const LifeV::Real& x, const LifeV::Real& y,
                               const LifeV::Real& z, const LifeV::ID& /*i*/ );

};

} // end namespace PyOrbLifeV


#endif /* _CUSTOM_FUNCTOR_HPP_ */
