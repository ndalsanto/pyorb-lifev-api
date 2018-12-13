/*
@author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
@date 27-11-2018
*/

#include "init_lifev.hpp"
#include "customFunctor.hpp"
#include <lifev/core/LifeV.hpp>
#include <mpi.h>
#include <Teuchos_RCP.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

// Cube's walls identifiers
const int BACK   = 101;
const int FRONT  = 102;
const int LEFT   = 103;
const int RIGHT  = 104;
const int BOTTOM = 105;
const int TOP    = 106;

namespace PyOrbLifeV
{

LifeV::Real zero_function( const LifeV::Real& /*t*/, const LifeV::Real& /*x*/, const LifeV::Real& /*y*/,
                const LifeV::Real& /*z*/, const LifeV::ID& /*i*/)
{
    return 0.;
}

LifeV::Real sourceFunction( const LifeV::Real& /*t*/, const LifeV::Real& /*x*/, const LifeV::Real& /*y*/,
                const LifeV::Real& /*z*/, const LifeV::ID& /*i*/)
{
    return 1.;
}

LifeV::Real mu1 = 1.0;
LifeV::Real mu2 = 1.0;
LifeV::Real mu3 = 1.0;

LifeV::Real diffusionFunction( const LifeV::Real& /*t*/, const LifeV::Real& x, const LifeV::Real& y,
                       const LifeV::Real& z, const LifeV::ID& /*i*/ )
{
    if( z < 0.5 && y < 0.5 ) return mu1;
    if( z > 0.5 && y < 0.5 ) return mu2;
    if( z < 0.5 && y > 0.5 ) return mu3;

    if( ( mu1 == 1.0 && mu2 == 0.0 && mu3 == 0.0 ) || ( mu1 == 0.0 && mu2 == 1.0 && mu3 == 0.0 ) || ( mu1 == 0.0 && mu2 == 0.0 && mu3 == 1.0 ) )
        return 0.0;

    return 1.;
}

LifeVSimulator::
LifeVSimulator( )
{
}

int
LifeVSimulator::
initialize( FemSpecifics& _femSpecifics )
{
    MPI_Comm * external_mpi_communicator = ( MPI_Comm * ) _femSpecifics.external_communicator;

#ifdef HAVE_MPI
    M_comm.reset( new Epetra_MpiComm ( *( external_mpi_communicator ) ) );
#else
    M_comm.reset( new Epetra_SerialComm );
#endif

    M_verbose = M_comm->MyPID( );

    if( M_verbose )
        std::cout << "LifeVSimulator::initialize M_dataFile" << std::endl;

    std::string dataFileName( _femSpecifics.datafile_path );
    M_dataFile.reset( new GetPot( dataFileName ) );

    initialize_simulation( );
}

int
LifeVSimulator::
finalize()
{
}

int
LifeVSimulator::
get_fem_dimension( )
{
    LifeV::VectorEpetra f( M_ETuFESpace->map(), LifeV::Unique );
    int my_size = f.epetraVector().MyLength( );

    return my_size;
}

int
LifeVSimulator::
initialize_simulation( )
{
    typedef LifeV::RegionMesh< LifeV::LinearTetra >                 mesh_Type;
    typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra >           FESpace_Type;
    typedef LifeV::ETFESpace< mesh_Type, LifeV::MapEpetra, 3, 1 >   uSpaceETA_Type;

    if( M_verbose )
        std::cout << "LifeVSimulator::initialize_simulation reading mesh" << std::endl;

    std::shared_ptr< mesh_Type > fullMeshPtr ( new mesh_Type ( M_comm ) );

    LifeV::MeshData meshData;
    meshData.setup ( *M_dataFile, "mesh");
    readMesh (*fullMeshPtr, meshData);

    LifeV::MeshPartitioner< mesh_Type > meshPart;

    meshPart.doPartition ( fullMeshPtr, M_comm );
    M_localMeshPtr = meshPart.meshPartition();

    // Clearing global mesh
    fullMeshPtr.reset();

    if( M_verbose )
        std::cout << "LifeVSimulator::initialize_simulation creating fe spaces" << std::endl;

    // Defining finite elements standard and ET spaces
    M_uFESpace.reset( new FESpace_Type ( M_localMeshPtr, (*M_dataFile)( "finite_element/degree", "P1" ), 1, M_comm ) );
    M_ETuFESpace.reset( new uSpaceETA_Type ( M_localMeshPtr, & ( M_uFESpace->refFE() ), & ( M_uFESpace->fe().geoMap() ), M_comm ) );

    M_A.reset( new LifeV::MatrixEpetra< LifeV::Real >( M_ETuFESpace->map(), 100 ) );
    M_f.reset( new LifeV::VectorEpetra( M_ETuFESpace->map(), LifeV::Unique ) );
    M_A->zero( );
    M_f->zero( );
}

int
LifeVSimulator::
build_stiffness_matrix( double * _param )
{
    mu1 = _param[0];
    mu2 = _param[1];
    mu3 = _param[2];

    std::cout << "Building stiffness matrix with parameters " << mu1 << " " << mu2 << " " << mu3 << std::endl;

    M_A.reset( new LifeV::MatrixEpetra< LifeV::Real >( M_ETuFESpace->map(), 100 ) );

    M_A->zero( );

    {
        using namespace std::placeholders;

        std::shared_ptr<customFunctor< LifeV::Real > >
                    diffusionFunctor ( new customFunctor< LifeV::Real >( diffusionFunction ) );


        {
            using namespace LifeV;
            using namespace ExpressionAssembly;

            integrate ( elements ( M_localMeshPtr ),
                    M_uFESpace->qr(),
                    M_ETuFESpace,
                    M_ETuFESpace,
                    value(1.0) * eval( diffusionFunctor, X ) * dot( grad( phi_j ), grad( phi_i ) )
                    )
                    >> M_A;
        }
    }

    M_A->globalAssemble( );

    applyBC( );
}

int
LifeVSimulator::
build_fem_vector( double * _param )
{
    mu1 = _param[0];
    mu2 = _param[1];
    mu3 = _param[2];

    std::shared_ptr< LifeV::VectorEpetra > f_repeated;
    f_repeated.reset( new LifeV::VectorEpetra( M_ETuFESpace->map(), LifeV::Repeated ) );
    f_repeated->zero();
    M_f.reset( new LifeV::VectorEpetra( M_ETuFESpace->map(), LifeV::Unique ) );
    M_f->zero();

    {
        using namespace LifeV;
        using namespace ExpressionAssembly;

        integrate (
                    elements ( M_localMeshPtr ),
                    M_uFESpace->qr(),
                    M_ETuFESpace,
                    value(100.0) * phi_i
                )
                >> f_repeated;
    }

    // f_repeated->globalAssemble( );
    *M_f += *f_repeated;

    applyBC( );
}

int
LifeVSimulator::
applyBC( )
{
    LifeV::BCHandler bcHandler;
    LifeV::BCFunctionBase ZeroBC ( zero_function );

    bcHandler.addBC( "Back",   BACK,   LifeV::Essential, LifeV::Full, ZeroBC, 1 );
    bcHandler.addBC( "Left",   LEFT,   LifeV::Essential, LifeV::Full, ZeroBC, 1 );
    bcHandler.addBC( "Top",    TOP,    LifeV::Essential, LifeV::Full, ZeroBC, 1 );

    bcHandler.addBC( "Front",  FRONT,  LifeV::Essential, LifeV::Full, ZeroBC, 1 );
    bcHandler.addBC( "Right",  RIGHT,  LifeV::Essential, LifeV::Full, ZeroBC, 1 );
    // bcHandler.addBC( "Bottom", BOTTOM, LifeV::Essential, LifeV::Full, ZeroBC, 1 );

    bcHandler.bcUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    LifeV::bcManage ( *M_A, *M_f, *M_uFESpace->mesh(), M_uFESpace->dof(), bcHandler, M_uFESpace->feBd(), 1.0, 0.0 );

    M_f->globalAssemble( );
    M_A->globalAssemble( );
}

int
LifeVSimulator::
perform_simulation( double * _param )
{
    build_stiffness_matrix( _param );
    build_fem_vector( _param );

    M_u.reset( new LifeV::VectorEpetra ( M_uFESpace->map(), LifeV::Unique ) );
    M_u->zero();

    LifeV::LinearSolver linearSolver ( M_comm );

    linearSolver.setOperator( M_A );

    Teuchos::RCP< Teuchos::ParameterList > aztecList = Teuchos::rcp ( new Teuchos::ParameterList );
    aztecList = Teuchos::getParametersFromXmlFile ( (*M_dataFile)( "solver/param_list_location", "SolverParamList.xml") );

    linearSolver.setParameters ( *aztecList );

    typedef LifeV::PreconditionerML                                         precML_type;
    typedef std::shared_ptr<precML_type>                                    precMLPtr_type;
    std::cout << "\t setting preconditioner ...";
    precML_type * precRawPtr;
    precRawPtr = new precML_type;
    // we set to look for the "fake" precMLL entry in order to set the default parameters
    // of ML preconditioner
    precRawPtr->setDataFromGetPot ( *M_dataFile, "precMLL" );
    std::shared_ptr<LifeV::Preconditioner> precPtr;
    precPtr.reset ( precRawPtr );
    std::cout << "done" << std::endl;

    linearSolver.setPreconditioner ( precPtr );
    linearSolver.setRightHandSide( M_f );
    linearSolver.solve( M_u );
}

double
LifeVSimulator::
compute_residual( )
{
    std::shared_ptr<LifeV::VectorEpetra> residual;
    residual.reset( new LifeV::VectorEpetra( M_ETuFESpace->map(), LifeV::Repeated ) );
    residual->zero();

    M_A->multiply( false, *M_u, *residual );
    *residual *= -1.;
    *residual += *M_f;

    return residual->norm2( );
}


}
