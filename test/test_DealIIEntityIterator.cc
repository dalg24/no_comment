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
    dealii_mesh->refine_global(2);

    Teuchos::RCP<DataTransferKit::DealIIAdjacencies<dim,spacedim>> adjacencies =
        Teuchos::rcp(new DataTransferKit::DealIIAdjacencies<dim,spacedim>(dealii_mesh));

    // Create a dtk entity iterator
    auto selectAll = [](DataTransferKit::Entity const &){return true;};
    DataTransferKit::EntityIterator dtk_entity_iterator =
        DataTransferKit::DealIIEntityIterator<dim,dim,spacedim>(
            dealii_mesh->begin_active(),
            dealii_mesh->begin_active(),
            dealii_mesh->end(),
            adjacencies.ptr(),
            selectAll
        );

    auto const iterator_begin = dtk_entity_iterator.begin();
    auto const iterator_end   = dtk_entity_iterator.end  ();
    BOOST_CHECK( dtk_entity_iterator == iterator_begin );
    BOOST_CHECK( dtk_entity_iterator != iterator_end   );
    std::set<DataTransferKit::EntityId> entity_ids;
    for (dtk_entity_iterator = iterator_begin;
        dtk_entity_iterator != iterator_end; ++dtk_entity_iterator)
    {
        entity_ids.emplace(dtk_entity_iterator->id());
        BOOST_CHECK_EQUAL( dtk_entity_iterator->ownerRank(), world.rank() );
    }
    BOOST_CHECK( dtk_entity_iterator == iterator_end   );
    BOOST_CHECK( dtk_entity_iterator != iterator_begin );

    // Check iterator size
    BOOST_CHECK_EQUAL( dtk_entity_iterator.size(), entity_ids.size() );

    // Print rank and set of entity ids owned
    std::cout<<world.rank()<<" : ";
    for (auto id : entity_ids)
        std::cout<<id<<"  ";
    std::cout<<"\n";

}
