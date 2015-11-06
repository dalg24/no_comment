#define BOOST_TEST_MODULE DealIIEntityShapeFunction
#include "main_included.cc"
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/format.hpp>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/fe_q.h>
#include <no_comment/DTK_DealIIEntitySet.h>
#include <no_comment/DTK_DealIIEntityShapeFunction.h>
#include <functional>

BOOST_AUTO_TEST_CASE( test_DealIIEntityShapeFunction )
{
    // Probably want to call templated function
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh with a single element
    boost::mpi::communicator world;
    Teuchos::RCP<DataTransferKit::DealIIMesh<dim,spacedim>> dealii_mesh =
        Teuchos::rcp(new DataTransferKit::DealIIMesh<dim,spacedim>(world));
    dealii::GridGenerator::hyper_rectangle(*dealii_mesh,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);

    // We use it to build a dtk entity on top of the element
    // adjacencies is a poor choice as for the name
    // it will change soon
    Teuchos::RCP<DataTransferKit::DealIIAdjacencies<dim,spacedim>> adjacencies =
        Teuchos::rcp(new DataTransferKit::DealIIAdjacencies<dim,spacedim>(dealii_mesh));

    // Distribute degrees of freedom
    Teuchos::RCP<dealii::DoFHandler<dim,spacedim>> dealii_dof_handler =
        Teuchos::rcp(new dealii::DoFHandler<dim,spacedim>(*dealii_mesh));
    Teuchos::RCP<dealii::FiniteElement<dim,spacedim>> dealii_fe =
        Teuchos::rcp(new dealii::FE_Q<dim>(1) );
    dealii_dof_handler->distribute_dofs(*dealii_fe);

    // Make a dtk entity shape function
    Teuchos::RCP<DataTransferKit::EntityShapeFunction> dtk_entity_shape_function =
        Teuchos::rcp(new DataTransferKit::DealIIEntityShapeFunction<dim,spacedim>(dealii_dof_handler));

    // Advance iterator to an element that is locally owned
    auto tria_iterator = dealii_mesh->begin_active();
    while (tria_iterator != dealii_mesh->end())
    {
        if (tria_iterator->is_locally_owned())
            break;
        ++tria_iterator;
    }

    // Ensure the iterator does not point to the end before dereferencing
    bool iterator_points_to_the_end = false;
    if ( tria_iterator != dealii_mesh->end() )
    {
        auto tria_accessor = *tria_iterator;

        // Create a dtk entity for the single volume element in the mesh
        DataTransferKit::Entity dtk_entity =
            DataTransferKit::DealIIEntity<dim,dim,spacedim>(tria_accessor, adjacencies);

        // Get support ids of the element
        Teuchos::Array<DataTransferKit::EntityId> dtk_entity_support_ids;
        dtk_entity_shape_function->entitySupportIds(
            dtk_entity, dtk_entity_support_ids );
        std::vector<dealii::types::global_dof_index> const dealii_dof_indices =
            { 0, 1, 2, 3, 4, 5, 6, 7, };
        BOOST_CHECK_EQUAL_COLLECTIONS(
            dtk_entity_support_ids.begin(), dtk_entity_support_ids.end(),
            dealii_dof_indices.cbegin(), dealii_dof_indices.cend() );

        double const tolerance = 1.0e-16;

        // Support nodes coordinates in the frame of reference and dof index
        std::vector<std::pair<std::vector<double>,int>> const support_nodes_xyz_dof =
            { { { 0.0, 0.0, 0.0 }, 0 },
              { { 1.0, 0.0, 0.0 }, 1 },
              { { 0.0, 1.0, 0.0 }, 2 },
              { { 1.0, 1.0, 0.0 }, 3 },
              { { 0.0, 0.0, 1.0 }, 4 },
              { { 1.0, 0.0, 1.0 }, 5 },
              { { 0.0, 1.0, 1.0 }, 6 },
              { { 1.0, 1.0, 1.0 }, 7 },
            };
        for (auto const & node : support_nodes_xyz_dof)
        {
            Teuchos::Array<double> values;
            dtk_entity_shape_function->evaluateValue(
                dtk_entity, node.first, values );
            int const dofs_per_cell = values.size();
            for (int dof = 0; dof < dofs_per_cell; ++dof)
                if (dof == node.second)
                    BOOST_CHECK_CLOSE( values[dof], 1.0, tolerance );
                else
                    BOOST_CHECK_CLOSE( values[dof], 0.0, tolerance );

            // I am too lazy to check all of them
            // gradients at the center of the cell will do
            Teuchos::Array<Teuchos::Array<double>> gradients;
            dtk_entity_shape_function->evaluateGradient(
                dtk_entity, node.first, gradients );
        }

        // Center of the element in the reference frame
        std::vector<double> center = { 0.5, 0.5, 0.5 };
        Teuchos::Array<double> values;
        dtk_entity_shape_function->evaluateValue(
            dtk_entity, center, values );
        for (auto v : values)
            BOOST_CHECK_CLOSE( v, 0.125, tolerance );

        Teuchos::Array<Teuchos::Array<double>> gradients;
        dtk_entity_shape_function->evaluateGradient(
            dtk_entity, center, gradients );
        // Here we are...
        std::vector<std::vector<double>> gradients_at_center =
             { {-0.25, -0.25, -0.25},
               { 0.25, -0.25, -0.25},
               {-0.25,  0.25, -0.25},
               { 0.25,  0.25, -0.25},
               {-0.25, -0.25,  0.25},
               { 0.25, -0.25,  0.25},
               {-0.25,  0.25,  0.25},
               { 0.25,  0.25,  0.25},
             };
        // Scale it and check
        // NB: has to match dimensions in the grid generator
        std::vector<double> scaling_factor = { 1.0, 2.0, 3.0 };
        for (int d = 0; d < 3; ++d)
            for (int i = 0; i < 8; ++i)
               BOOST_CHECK_CLOSE(
                  gradients_at_center[i][d] / scaling_factor[d],
                  gradients[i][d], tolerance );
    } else {
        iterator_points_to_the_end = true;
    }
    // Ensure at least one process went through the if statement
    BOOST_CHECK( !boost::mpi::all_reduce(world, iterator_points_to_the_end, std::logical_and<bool>()) );
    
}
