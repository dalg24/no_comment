#define BOOST_TEST_MODULE DealIIEntity
#define BOOST_TEST_MAIN
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <boost/test/unit_test.hpp>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/mpi.h>
#include <no_comment/DTK_DealIIEntity.h>

#include "MPIFixture.cc"

BOOST_FIXTURE_TEST_CASE( test_DealIIEntity_cell, MPIFixture )
{
    // Probably want to call templated function
    int const structdim = 3;
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh
    DataTransferKit::DealIIMesh<dim,spacedim> dealii_mesh(world);
    dealii::GridGenerator::hyper_rectangle(dealii_mesh,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);

    auto tria_iterator = dealii_mesh.begin_active();
    while (tria_iterator != dealii_mesh.end())
    {
        if (tria_iterator->is_locally_owned()) break;
        ++tria_iterator;
    }

    bool iterator_points_to_the_end = false;
    if ( tria_iterator != dealii_mesh.end() )
    {
        auto tria_accessor = *tria_iterator;

        // Create a dtk entity for the single volume element in the mesh
        DataTransferKit::Entity dtk_entity =
            DataTransferKit::DealIIEntity<structdim,dim,spacedim>(tria_accessor);

        // Print out the entity.
        Teuchos::RCP<Teuchos::FancyOStream>
            fancy_os = Teuchos::VerboseObjectBase::getDefaultOStream();
        dtk_entity.describe( *fancy_os );

        // Check topological and physical dimensions
        BOOST_CHECK_EQUAL( dtk_entity.topologicalDimension(), dim      );
        BOOST_CHECK_EQUAL( dtk_entity.physicalDimension()   , spacedim );

        // Check owner rank
        BOOST_CHECK_EQUAL( dtk_entity.ownerRank(), world.rank() );

        // Bounding box
        Teuchos::Tuple<double,6> bounding_box;
        dtk_entity.boundingBox(bounding_box);

        double const tolerance = 1.0e-15;
        BOOST_CHECK_SMALL( -1.0 - bounding_box[0], tolerance );
        BOOST_CHECK_SMALL( -2.0 - bounding_box[1], tolerance );
        BOOST_CHECK_SMALL( -3.0 - bounding_box[2], tolerance );
        BOOST_CHECK_SMALL(  0.0 - bounding_box[3], tolerance );
        BOOST_CHECK_SMALL(  0.0 - bounding_box[4], tolerance );
        BOOST_CHECK_SMALL(  0.0 - bounding_box[5], tolerance );

        // Block id and boundary id
        BOOST_CHECK(  dtk_entity.inBlock(0) );
        BOOST_CHECK( !dtk_entity.inBlock(1) );
        BOOST_CHECK( !dtk_entity.inBlock(2) );

        BOOST_CHECK(  dtk_entity.onBoundary(0) );
        BOOST_CHECK(  dtk_entity.onBoundary(1) );
        BOOST_CHECK(  dtk_entity.onBoundary(2) );
        BOOST_CHECK(  dtk_entity.onBoundary(3) );
        BOOST_CHECK(  dtk_entity.onBoundary(4) );
        BOOST_CHECK(  dtk_entity.onBoundary(5) );
        BOOST_CHECK( !dtk_entity.onBoundary(6) );
    } else {
        iterator_points_to_the_end = true;
    }
    // Ensure at least one process went through the if statement
    BOOST_CHECK( !boost::mpi::all_reduce(world, iterator_points_to_the_end, std::logical_and<bool>()) );
}

// TODO: FIX THE NODE TEST
BOOST_FIXTURE_TEST_CASE( test_DealIIEntity_node, MPIFixture )
{
    // Probably want to call templated function
    int const structdim = 0;
    int const dim = 3;
    int const spacedim = 3;

    // Build a distributed mesh
    DataTransferKit::DealIIMesh<dim,spacedim> dealii_mesh(world);
    dealii::GridGenerator::hyper_rectangle(dealii_mesh,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);

    // Advance iterator to a cell that is locally owned
    auto tria_iterator = dealii_mesh.begin_active();
    while (tria_iterator != dealii_mesh.end())
    {
        if (tria_iterator->is_locally_owned()) break;
        ++tria_iterator;
    }

    auto tria_accessor = *tria_iterator;

    auto vertex_to_cell = dealii::GridTools::vertex_to_cell_map(dealii_mesh);
    auto local_to_global_vertex_id =
      dealii::GridTools::compute_local_to_global_vertex_index_map(dealii_mesh);
    Teuchos::RCP<std::vector<std::set<
      typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>>
      teuchos_vertex_to_cell = Teuchos::rcpFromRef(vertex_to_cell);
    Teuchos::RCP<std::map<unsigned int, unsigned long long int>>
      teuchos_local_to_global_vertex_id = Teuchos::rcpFromRef(local_to_global_vertex_id);

    bool boundary[8][6] =
    {
      {true,  false, true,  false, true,  false},
      {false, true,  true,  false, true,  false},
      {true,  false, false, true,  true,  false},
      {false, true,  false, true,  true,  false},
      {true,  false, true,  false, false, true },
      {false, true,  true,  false, false, true },
      {true,  false, false, true,  false, true },
      {false, true,  false, true,  false, true },
    };

    // Create a dtk entity for some vertices
    for (unsigned int i=0; i<8; ++i)
    {
      auto vertex_iterator = tria_accessor.vertex_iterator(i);
      DataTransferKit::Entity dtk_entity =
        DataTransferKit::DealIIEntity<structdim,dim,spacedim>(*vertex_iterator, teuchos_vertex_to_cell,
            teuchos_local_to_global_vertex_id);
      // Print out the entity.
      Teuchos::RCP<Teuchos::FancyOStream>
        fancy_os = Teuchos::VerboseObjectBase::getDefaultOStream();
      dtk_entity.describe( *fancy_os );

      // Check topological and physical dimensions
      BOOST_CHECK_EQUAL( dtk_entity.topologicalDimension(), dim      );
      BOOST_CHECK_EQUAL( dtk_entity.physicalDimension()   , spacedim );

      // Bounding box
      Teuchos::Tuple<double,6> bounding_box;
      dtk_entity.boundingBox(bounding_box);

      double const tolerance = 1.0e-15;
      BOOST_CHECK_SMALL( vertex_iterator->center()[0] - bounding_box[0], tolerance );
      BOOST_CHECK_SMALL( vertex_iterator->center()[1] - bounding_box[1], tolerance );
      BOOST_CHECK_SMALL( vertex_iterator->center()[2] - bounding_box[2], tolerance );
      BOOST_CHECK_SMALL( vertex_iterator->center()[0] - bounding_box[3], tolerance );
      BOOST_CHECK_SMALL( vertex_iterator->center()[1] - bounding_box[4], tolerance );
      BOOST_CHECK_SMALL( vertex_iterator->center()[2] - bounding_box[5], tolerance );

      // Block id and boundary id
      BOOST_CHECK(  dtk_entity.inBlock(0) );
      BOOST_CHECK( !dtk_entity.inBlock(1) );
      BOOST_CHECK( !dtk_entity.inBlock(2) );

      for (unsigned int j=0; j<6; ++j)
      BOOST_CHECK(  dtk_entity.onBoundary(j) == boundary[i][j] );
    }
}
