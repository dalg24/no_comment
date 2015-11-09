#define BOOST_TEST_MODULE DealIIEntityIterator
#include "main_included.cc"
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

BOOST_AUTO_TEST_CASE( test_DealIIEntitySet )
{
    // Probably want to call templated function
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh
    boost::mpi::communicator world;
    Teuchos::RCP<DataTransferKit::DealIIMesh<dim,spacedim>> dealii_mesh =
        Teuchos::rcp(new DataTransferKit::DealIIMesh<dim,spacedim>(world));

    dealii::GridGenerator::hyper_rectangle(*dealii_mesh,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);

    Teuchos::RCP<DataTransferKit::DealIIAdjacencies<dim,spacedim>> adjacencies =
        Teuchos::rcp(new DataTransferKit::DealIIAdjacencies<dim,spacedim>(dealii_mesh));

    // Create a dtk entity set
    auto selectAll = [](DataTransferKit::Entity const &){return true;};
    DataTransferKit::EntityIterator dtk_entity_iterator =
        DataTransferKit::DealIIEntityIterator<dim,dim,spacedim>(
            dealii_mesh->begin_active(),
            dealii_mesh->begin_active(),
            dealii_mesh->end(),
            adjacencies.ptr(),
            selectAll
        );
}
