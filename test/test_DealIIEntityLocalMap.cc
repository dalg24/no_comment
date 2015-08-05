#define BOOST_TEST_MODULE DealIIEntityLocalMap
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <no_comment/DTK_DealIIEntity.h>
#include <no_comment/DTK_DealIIEntityLocalMap.h>

BOOST_AUTO_TEST_CASE( test_DealIIEntityLocalMap )
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
        std::make_shared<dealii::TriaAccessor<structdim,dim,spacedim>>(
            *tria.begin_active() );

    // Create a dtk entity for the single volume element in the mesh
    DataTransferKit::Entity dtk_entity =
        DealIIEntity<structdim,dim,spacedim>(tria_iterator);

}

