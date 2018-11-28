/*
@author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
@date 27-11-2018
*/

#ifndef __INIT_LIFEV_PYORB_API__
#define __INIT_LIFEV_PYORB_API__

#include <iostream>
#include <memory>
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>

namespace PyOrbLifeV
{

class FemSpecifics
{
public:
      char * model;
      char * datafile_path;
      MPI_Comm * external_communicator;
      double * u;
      double * A;
      double * f;
};

class LifeVSimulator
{

public:

    LifeVSimulator( );

    ~LifeVSimulator( ) { };

    // methods to initialize/finalize LifeV
    int initialize( FemSpecifics& _femSpecifics );

    int finalize();

    // methods to run a LifeV simulation
    int initialize_simulation( );

    int perform_simulation( double * _param );

    int get_fem_dimension( );

    int build_stiffness_matrix( double * _param );

    int build_fem_vector( double * _param );

    std::shared_ptr< LifeV::VectorEpetra > get_solution( )
    {
        return M_u;
    }

    // members to initialize LifeV
    std::shared_ptr< Epetra_Comm >                      M_comm;

    // members to initialize a LifeV simulation
    std::shared_ptr<GetPot>                                           M_dataFile;
    std::shared_ptr< LifeV::RegionMesh< LifeV::LinearTetra > >        M_localMeshPtr;
    std::shared_ptr< LifeV::FESpace< LifeV::RegionMesh< LifeV::LinearTetra >, LifeV::MapEpetra > >            M_uFESpace;
    std::shared_ptr< LifeV::ETFESpace< LifeV::RegionMesh< LifeV::LinearTetra >, LifeV::MapEpetra, 3, 1 > >    M_ETuFESpace;
    std::shared_ptr< LifeV::VectorEpetra > M_u;
    std::shared_ptr< LifeV::VectorEpetra > M_f;
    std::shared_ptr< LifeV::MatrixEpetra< LifeV::Real > > M_A;
};

}



#endif
