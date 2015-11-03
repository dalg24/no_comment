#define BOOST_TEST_MODULE DealIIEntityIterator
#define BOOST_TEST_MAIN
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/format.hpp>
#include <boost/mpi.hpp>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/mpi.h>
#include <no_comment/DTK_DealIIEntityIterator.h>
#include <functional>

#include "MPIFixture.cc"

BOOST_FIXTURE_TEST_CASE( test_DealIIEntitySet, MPIFixture )
{
    // Probably want to call templated function
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh
    Teuchos::RCP<DataTransferKit::DealIIMesh<dim,spacedim>> dealii_mesh =
        Teuchos::rcp(new DataTransferKit::DealIIMesh<dim,spacedim>(world));

    dealii::GridGenerator::hyper_rectangle(*dealii_mesh,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);

    // Create a dtk entity set
    auto select_all = [](DataTransferKit::Entity const &){return true;};
//    DataTransferKit::EntityIterator dtk_entity_iterator =
//        DealIIEntityIterator<dim,spacedim>(
//            dealii_mesh->begin_active(),
//            dealii_mesh->begin_active(),
//            dealii_mesh->end(),
//            Teuchos::Ptr<dealii::parallel::distributed::Triangulation<dim,spacedim>>(dealii_mesh.get()),
//            Teuchos::Ptr<DealIIAdjacencies<dim,spacedim>>(),
//            select_all
//        );
}
