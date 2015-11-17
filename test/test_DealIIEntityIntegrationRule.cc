#define BOOST_TEST_MODULE DealIIEntityIntegrationRule
#include "main_included.cc"
#include <boost/test/unit_test.hpp>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <no_comment/DTK_DealIIEntity.h>
#include <no_comment/DTK_DealIIEntityIntegrationRule.h>

BOOST_AUTO_TEST_CASE( test_DealIIEntityIntegrationRule )
{
    // Probably want to call templated function
    int const structdim = 3;
    int const dim = 3;
    int const spacedim = 3;

    boost::mpi::communicator world;

    // Build a mesh
    Teuchos::RCP<DataTransferKit::DealIIMesh<dim,spacedim>> dealii_mesh =
        Teuchos::rcp(new DataTransferKit::DealIIMesh<dim,spacedim>(world));
    dealii::GridGenerator::hyper_rectangle(*dealii_mesh,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);

    Teuchos::RCP<DataTransferKit::DealIIAdjacencies<dim,spacedim>> adjacencies =
        Teuchos::rcp(new DataTransferKit::DealIIAdjacencies<dim,spacedim>(dealii_mesh));

    auto dealii_tria_iterator =
            dealii_mesh->begin_active();

    // Create a dtk entity for the single volume element in the mesh
    DataTransferKit::Entity dtk_entity =
        DataTransferKit::DealIIEntity<structdim,dim,spacedim>(*dealii_tria_iterator, adjacencies.ptr());

    // Create an integration rule
    Teuchos::RCP<DataTransferKit::EntityIntegrationRule> dtk_entity_integration_rule =
        Teuchos::rcp(new DataTransferKit::DealIIEntityIntegrationRule());

    Teuchos::Array<Teuchos::Array<double>> points;
    Teuchos::Array<double> weights;
    for (int order = 0; order < 3; ++ order) 
    {
        std::cout<<order<<"\n";
        dtk_entity_integration_rule->getIntegrationRule(
            dtk_entity,
            order,
            points,
            weights );
        size_t const n = points.size();
        BOOST_CHECK_EQUAL( n, weights.size() );
        for (size_t q = 0; q < n; ++q)
        {
            BOOST_CHECK_EQUAL( points[q].size(), 3 );
            for (auto const & xyz : points[q])
                std::cout<<"  "<<xyz;
            std::cout<<" ("<<weights[q]<<")\n";

        }
    }

}

