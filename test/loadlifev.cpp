#include <iostream>
#include "test_function.hpp"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChronoManager.hpp>

int main (int argc, char **argv)
{
//	std::cout << "LIFEV version: " << LIFEV_VERSION_MAJOR
//		  << "." << LIFEV_VERSION_MINOR << std::endl;

    std::cout << "Hello world!" << std::endl;

    LifeV::LifeChrono globalChrono;


    test_function( 3 );

    return 0;
}
