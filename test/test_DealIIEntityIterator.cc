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
            adjacencies.ptr(),
            selectAll
        );

    auto const iterator_begin = dtk_entity_iterator.begin();
    auto const iterator_end   = dtk_entity_iterator.end  ();
    BOOST_CHECK( dtk_entity_iterator == iterator_begin );
    BOOST_CHECK( dtk_entity_iterator != iterator_end   );
    std::map<DataTransferKit::EntityId,int> entities;
    for (dtk_entity_iterator = iterator_begin;
        dtk_entity_iterator != iterator_end; ++dtk_entity_iterator)
    {
        entities.emplace(
            std::pair<DataTransferKit::EntityId,int>(
                dtk_entity_iterator->id(),dtk_entity_iterator->ownerRank() )
        );
//        BOOST_CHECK_EQUAL( dtk_entity_iterator->ownerRank(), world.rank() );
    }
    BOOST_CHECK( dtk_entity_iterator == iterator_end   );
    BOOST_CHECK( dtk_entity_iterator != iterator_begin );

    // Check iterator size
    BOOST_CHECK_EQUAL( dtk_entity_iterator.size(), entities.size() );

    // Print rank(number of entities) : set of entity id(owner rank)
    std::cout<<world.rank()<<"("<<entities.size()<<") : ";
    for (auto entity : entities)
        std::cout<<entity.first<<"("<<entity.second<<")  ";
    std::cout<<"\n";


    // Select locally owned elems
    int const rank = world.rank();
    auto selectLocallyOwned =
        [rank](DataTransferKit::Entity const & entity)
        {
            return entity.ownerRank() == rank;
        };
     dtk_entity_iterator =
        DataTransferKit::DealIIEntityIterator<dim,dim,spacedim>(
            adjacencies.ptr(),
            selectLocallyOwned
        );
     size_t const global_size = boost::mpi::all_reduce(
         world, dtk_entity_iterator.size(), std::plus<size_t>() );
     BOOST_CHECK_EQUAL( global_size, 64 );

}
