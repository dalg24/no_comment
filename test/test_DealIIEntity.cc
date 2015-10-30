#define BOOST_TEST_MODULE DealIIEntity
#define BOOST_TEST_MAIN
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <boost/test/unit_test.hpp>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/mpi.h>
#include <no_comment/DTK_DealIIEntity.h>

#include "MPIFixture.cc"
BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE( test_DealIIEntity_cell )
{
    // Probably want to call templated function
    int const structdim = 3;
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh
    dealii::parallel::distributed::Triangulation<dim,spacedim> tria(MPI_COMM_WORLD);
    dealii::GridGenerator::hyper_rectangle(tria,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);

    auto tria_iterator =
            *tria.begin_active();

    // Create a dtk entity for the single volume element in the mesh
    DataTransferKit::Entity dtk_entity =
        DealIIEntity<structdim,dim,spacedim>(tria_iterator);

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
}

BOOST_AUTO_TEST_CASE( test_DealIIEntity_node )
{
  // Probably want to call templated function
  int const structdim = 0;
  int const dim = 3;
  int const spacedim = 3;

  // Build a distributed mesh
  dealii::parallel::distributed::Triangulation<dim,spacedim> tria(MPI_COMM_WORLD);
  dealii::GridGenerator::hyper_rectangle(tria,
      dealii::Point<spacedim>(-1.0, -2.0, -3.0),
      dealii::Point<spacedim>( 0.0,  0.0,  0.0),
      true);

  auto tria_iterator = *tria.begin_active();
  auto vertex_to_cell = dealii::GridTools::vertex_to_cell_map(tria);
  auto local_to_global_vertex_id = 
    dealii::GridTools::compute_local_to_global_vertex_index_map(tria);
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
    auto vertex_iterator = tria_iterator.vertex_iterator(i);
    DataTransferKit::Entity dtk_entity =
      DealIIEntity<structdim,dim,spacedim>(*vertex_iterator, teuchos_vertex_to_cell, 
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
