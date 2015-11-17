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
    int const n = 1;
    dealii_mesh->refine_global(n);
    int const elems = 8;
    int const nodes_per_boundary = 9;

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
    int const rank = world.rank();
    DataTransferKit::PredicateFunction localyOwned =
        [rank](DataTransferKit::Entity const & entity)
        { return (entity.ownerRank() == rank); };
    DataTransferKit::EntityIterator dtk_elem_iterator =
        dtk_entity_set->entityIterator(dim, localyOwned);
    BOOST_CHECK_EQUAL(
        boost::mpi::all_reduce(world, dtk_elem_iterator.size(), std::plus<size_t>()),
        elems );

    // Get iterator over the surface nodes
    int const boundary_id = 0;
    DataTransferKit::PredicateFunction onBoundaryAndLocalyOwned =
        [rank,boundary_id](DataTransferKit::Entity const & entity)
        { return entity.onBoundary(boundary_id) && (entity.ownerRank() == rank); };
    DataTransferKit::EntityIterator dtk_node_iterator =
        dtk_entity_set->entityIterator(0, onBoundaryAndLocalyOwned);
    BOOST_CHECK_EQUAL(
        boost::mpi::all_reduce(world, dtk_node_iterator.size(), std::plus<size_t>()),
        nodes_per_boundary );

    // Check the bounding box
    double const tolerance = 1.0e-15;
    Teuchos::Tuple<double,6> bounding_box;
    dtk_entity_set->globalBoundingBox(bounding_box);
    BOOST_CHECK_SMALL( -1.0 - bounding_box[0], tolerance );
    BOOST_CHECK_SMALL( -2.0 - bounding_box[1], tolerance );
    BOOST_CHECK_SMALL( -3.0 - bounding_box[2], tolerance );
    BOOST_CHECK_SMALL(  0.0 - bounding_box[3], tolerance );
    BOOST_CHECK_SMALL(  0.0 - bounding_box[4], tolerance );
    BOOST_CHECK_SMALL(  0.0 - bounding_box[5], tolerance );

    // Get adjacencies
    Teuchos::Array<DataTransferKit::Entity> adjacent_elem;

}
