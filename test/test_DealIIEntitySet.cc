#define BOOST_TEST_MODULE DealIIEntitySet
#include "main_included.cc"
#include <no_comment/DTK_DealIIEntitySet.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <boost/test/unit_test.hpp>
#include <boost/mpi.hpp>

BOOST_AUTO_TEST_CASE( test_DealIIEntitySet )
{
    // Probably want to call templated function
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh
    boost::mpi::communicator world;
    Teuchos::RCP<DataTransferKit::DealIIMesh<dim,spacedim>> dealii_mesh =
        Teuchos::rcp(new DataTransferKit::DealIIMesh<dim,spacedim>(world));

    dealii::GridGenerator::hyper_rectangle(
        *dealii_mesh,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);

    // Create a dtk entity set
    Teuchos::RCP<DataTransferKit::EntitySet> dtk_entity_set =
        Teuchos::rcp(new DataTransferKit::DealIIEntitySet<dim,spacedim>(dealii_mesh));

    // Get the communicator
    auto comm = dtk_entity_set->communicator();
    BOOST_ASSERT( world.rank() == comm->getRank() );
    BOOST_ASSERT( world.size() == comm->getSize() );

    // Check physical dimension
    BOOST_ASSERT( dtk_entity_set->physicalDimension() == spacedim );

    // Get iterator over the volume elements
    DataTransferKit::EntityIterator dtk_entity_iterator =
        dtk_entity_set->entityIterator(dim);
    int const global_size =
        boost::mpi::all_reduce(world, dtk_entity_iterator.size(), std::plus<double>());
    BOOST_CHECK_EQUAL( global_size, 1 );

    Teuchos::Tuple<double,6> bounding_box;
    dtk_entity_set->globalBoundingBox(bounding_box);
    BOOST_CHECK_EQUAL( bounding_box[0], -1.0 );
    BOOST_CHECK_EQUAL( bounding_box[1], -2.0 );
    BOOST_CHECK_EQUAL( bounding_box[2], -3.0 );
    BOOST_CHECK_EQUAL( bounding_box[3],  0.0 );
    BOOST_CHECK_EQUAL( bounding_box[4],  0.0 );
    BOOST_CHECK_EQUAL( bounding_box[5],  0.0 );

    // Get adjacencies
    Teuchos::Array<DataTransferKit::Entity> adjacent_elem;

}
