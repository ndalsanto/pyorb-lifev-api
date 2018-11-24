#include <iostream>
#include "test_function.hpp"
// #include <rb/reduced_basis/fem/FemProblem.hpp>
        
int main (int argc, char **argv)
{
//	std::cout << "LIFEV version: " << LIFEV_VERSION_MAJOR
//		  << "." << LIFEV_VERSION_MINOR << std::endl;

	std::cout << "Ciao" << std::endl;

	test_function( 3 );

	return 0; 
}

