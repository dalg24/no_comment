#define BOOST_TEST_MODULE DealIIEntity
#define BOOST_TEST_MAIN
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <boost/test/unit_test.hpp>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <no_comment/DTK_DealIIEntity.h>

BOOST_AUTO_TEST_CASE( test_DealIIEntity )
{
    // Probably want to call templated function
    int const structdim = 3;
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh
    dealii::Triangulation<dim,spacedim> tria;
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

