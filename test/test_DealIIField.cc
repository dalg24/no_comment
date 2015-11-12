#define BOOST_TEST_MODULE DealIIField
#include "main_included.cc"
#include <no_comment/DTK_DealIIField.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <boost/test/unit_test.hpp>
#include <boost/mpi.hpp>

BOOST_AUTO_TEST_CASE( test_DealIIEntitySet )
{
    // Probably want to call templated function
    int const dim = 3;
    int const spacedim = 3;

    // Build a mesh
    boost::mpi::communicator world;
    Teuchos::RCP<DataTransferKit::DealIIMesh<dim,spacedim>> dealii_mesh =
        Teuchos::rcp(new DataTransferKit::DealIIMesh<dim,spacedim>(world));

    dealii::GridGenerator::hyper_rectangle(
        *dealii_mesh,
        dealii::Point<spacedim>(-1.0, -2.0, -3.0),
        dealii::Point<spacedim>( 0.0,  0.0,  0.0),
        true);
    dealii_mesh->refine_global(2);

//    // We use it to build a dtk entity on top of the element
//    // adjacencies is a poor choice as for the name
//    // it will change soon
//    Teuchos::RCP<DataTransferKit::DealIIAdjacencies<dim,spacedim>> adjacencies =
//        Teuchos::rcp(new DataTransferKit::DealIIAdjacencies<dim,spacedim>(dealii_mesh));

    // Distribute degrees of freedom
    Teuchos::RCP<dealii::DoFHandler<dim,spacedim>> dealii_dof_handler =
        Teuchos::rcp(new dealii::DoFHandler<dim,spacedim>(*dealii_mesh));
    Teuchos::RCP<dealii::FiniteElement<dim,spacedim>> dealii_fe =
        Teuchos::rcp(new dealii::FE_Q<dim>(1) );
    dealii_dof_handler->distribute_dofs(*dealii_fe);

    // Build a vector
    dealii::IndexSet const & locally_owned_dofs = dealii_dof_handler->locally_owned_dofs();
    dealii::IndexSet locally_relevant_dofs;
    dealii::DoFTools::extract_locally_relevant_dofs(*dealii_dof_handler, locally_relevant_dofs);
    Teuchos::RCP<DataTransferKit::DealIIVector> dealii_vector =
        Teuchos::rcp(new DataTransferKit::DealIIVector(locally_owned_dofs, locally_relevant_dofs, world));

    // Create a dtk field
    Teuchos::RCP<DataTransferKit::Field> dtk_field =
        Teuchos::rcp(new DataTransferKit::DealIIField(dealii_vector));

    // Check dimension of the field
    BOOST_CHECK_EQUAL( dtk_field->dimension(), 1 );

    // Get local support ids
    Teuchos::Array<DataTransferKit::SupportId> support_ids =
        dtk_field->getLocalSupportIds();
    BOOST_CHECK_EQUAL( support_ids.size(), dealii_vector->local_size() );

    // Print rank(number of ids) : support ids
    std::cout<<world.rank()<<" ("<<support_ids.size()<<") :";
    for (auto const & id : support_ids)
        std::cout<<id<<"  ";
    std::cout<<"\n";

    *dealii_vector = 3.14;



//    // Get the communicator
//    auto comm = dtk_entity_set->communicator();
//    BOOST_ASSERT( world.rank() == comm->getRank() );
//    BOOST_ASSERT( world.size() == comm->getSize() );
//
//    // Check physical dimension
//    BOOST_ASSERT( dtk_entity_set->physicalDimension() == spacedim );
//
//    // Get iterator over the volume elements
//    DataTransferKit::EntityIterator dtk_entity_iterator =
//        dtk_entity_set->entityIterator(dim);
//    int const global_size =
//        boost::mpi::all_reduce(world, dtk_entity_iterator.size(), std::plus<double>());
//    BOOST_CHECK_EQUAL( global_size, 1 );
//
//    Teuchos::Tuple<double,6> bounding_box;
//    dtk_entity_set->globalBoundingBox(bounding_box);
//    BOOST_CHECK_EQUAL( bounding_box[0], -1.0 );
//    BOOST_CHECK_EQUAL( bounding_box[1], -2.0 );
//    BOOST_CHECK_EQUAL( bounding_box[2], -3.0 );
//    BOOST_CHECK_EQUAL( bounding_box[3],  0.0 );
//    BOOST_CHECK_EQUAL( bounding_box[4],  0.0 );
//    BOOST_CHECK_EQUAL( bounding_box[5],  0.0 );
//
//    // Get adjacencies
//    Teuchos::Array<DataTransferKit::Entity> adjacent_elem;

}
