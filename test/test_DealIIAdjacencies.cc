#define BOOST_TEST_MODULE DealIIAdjacencies
#define BOOST_TEST_MAIN
#include <no_comment/DTK_DealIIAdjacencies.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/mpi.h>
#include <boost/test/unit_test.hpp>

#include "MPIFixture.cc"

BOOST_FIXTURE_TEST_CASE( test_DealIIAdjacencies, MPIFixture )
{
  // Probably want to call templated function
  int const dim = 3;
  int const spacedim = 3;

  // Build a mesh
  Teuchos::RCP<dealii::parallel::distributed::Triangulation<dim,spacedim>> tria =
    Teuchos::rcp(new dealii::parallel::distributed::Triangulation<dim,spacedim>(world));

  dealii::GridGenerator::hyper_rectangle(*tria,
      dealii::Point<spacedim>(-1.0, -2.0, -3.0),
      dealii::Point<spacedim>( 0.0,  0.0,  0.0),
      true);
  tria->refine_global(3);

  DealIIAdjacencies<dim,spacedim> dealii_adjacencies(tria);
}
