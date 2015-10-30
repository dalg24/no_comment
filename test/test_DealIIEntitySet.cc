#define BOOST_TEST_MODULE DealIIEntitySet
#define BOOST_TEST_MAIN
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/format.hpp>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/mpi.h>
#include <no_comment/DTK_DealIIEntitySet.h>

#include "MPIFixture.cc"
BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE( test_DealIIEntitySet )
{
    // Probably want to call templated function
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh
    Teuchos::RCP<dealii::parallel::distributed::Triangulation<dim,spacedim>> dealii_mesh =
        Teuchos::rcp(
            new dealii::parallel::distributed::Triangulation<dim,spacedim>(MPI_COMM_WORLD)
        );
    dealii::GridGenerator::hyper_rectangle(*dealii_mesh,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);

    // Create a dtk entity set
    Teuchos::RCP<DataTransferKit::EntitySet> dtk_entity_set =
        Teuchos::rcp(new DealIIEntitySet<dim,spacedim>(dealii_mesh));

    // Get the communicator
    auto comm = dtk_entity_set->communicator();
    std::cout<<boost::format("hello world from processor %d out of %d\n")
        % comm->getRank()
        % comm->getSize()
        ;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    BOOST_ASSERT( rank == comm->getRank() );
    BOOST_ASSERT( size == comm->getSize() );

    // Check physical dimension
    BOOST_ASSERT( dtk_entity_set->physicalDimension() == spacedim );

    // Get adjacencies
    Teuchos::Array<DataTransferKit::Entity> adjacent_elem;

}
