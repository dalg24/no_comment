#define BOOST_TEST_MODULE DealIIAdjacencies
#define BOOST_TEST_MAIN
#include <no_comment/DTK_DealIIAdjacencies.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test_DealIIAdjacencies )
{
    // Probably want to call templated function
    int const structdim = 3;
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh
    auto tria = std::make_shared<dealii::Triangulation<dim,spacedim>>();
    dealii::GridGenerator::hyper_rectangle(*tria,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);
    tria->refine_global(3);

    DealIIAdjacencies<dim,spacedim> dealii_adjacencies(tria);

}
