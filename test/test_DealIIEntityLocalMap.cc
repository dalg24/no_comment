#define BOOST_TEST_MODULE DealIIEntityLocalMap
#include "main_included.cc"
#include <boost/test/unit_test.hpp>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q1.h>
#include <no_comment/DTK_DealIIEntity.h>
#include <no_comment/DTK_DealIIEntityLocalMap.h>

BOOST_AUTO_TEST_CASE( test_DealIIEntityLocalMap_elem )
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

    // Create a local map
    auto dealii_mapping =
        std::make_shared<dealii::MappingQ1<dim,spacedim>>();
    Teuchos::RCP<DataTransferKit::EntityLocalMap> dtk_entity_local_map =
        Teuchos::rcp(new DataTransferKit::DealIIEntityLocalMap<structdim,dim,spacedim>(dealii_mapping) );

    // Measure
    double const percent_tolerance = 1.0e-10;
    BOOST_CHECK_CLOSE( dtk_entity_local_map->measure(dtk_entity), 6.0, percent_tolerance );

    // Centroid
    Teuchos::Array<double> centroid(spacedim, 0.0);
    dtk_entity_local_map->centroid(dtk_entity, centroid());
    BOOST_CHECK_CLOSE( centroid[0], -0.5, percent_tolerance );
    BOOST_CHECK_CLOSE( centroid[1], -1.0, percent_tolerance );
    BOOST_CHECK_CLOSE( centroid[2], -1.5, percent_tolerance );

    // Safe to map to reference frame
    Teuchos::Array<double> good_point(spacedim);
    good_point[0] = -0.15;
    good_point[1] = -0.1;
    good_point[2] = -1.2;
    Teuchos::Array<double> bad_point(spacedim);
    bad_point[0] = 0.3;
    bad_point[1] = 2.1;
    bad_point[2] = 0.1;
    BOOST_CHECK(  dtk_entity_local_map->isSafeToMapToReferenceFrame(dtk_entity, good_point()) ); 
    BOOST_CHECK( !dtk_entity_local_map->isSafeToMapToReferenceFrame(dtk_entity, bad_point ()) ); 

    // Map to reference frame
    Teuchos::Array<double> ref_good_point(spacedim);
    dtk_entity_local_map->mapToReferenceFrame(dtk_entity, good_point(), ref_good_point());
    for (unsigned int vertex = 0; vertex < dealii::GeometryInfo<dim>::vertices_per_cell; ++vertex)
        std::cout<<"vertex "<<vertex<<": "<<dealii_tria_iterator->vertex(vertex)<<"\n";
    // NOTE: I did ``1 - eta'' because of the cell orientation
    BOOST_CHECK_CLOSE( ref_good_point[0], 1.0 - 0.15, percent_tolerance );
    BOOST_CHECK_CLOSE( ref_good_point[1], 1.0 - 0.05, percent_tolerance );
    BOOST_CHECK_CLOSE( ref_good_point[2], 1.0 - 0.4 , percent_tolerance );

    // Check point inclusion 
    Teuchos::Array<double> ref_bad_point(spacedim);
    ref_bad_point[0] = -0.3;
    ref_bad_point[1] = 0.7;
    ref_bad_point[2] = 0.2;
    BOOST_CHECK(  dtk_entity_local_map->checkPointInclusion(dtk_entity, ref_good_point()) );
    BOOST_CHECK( !dtk_entity_local_map->checkPointInclusion(dtk_entity, ref_bad_point ()) );
    
    // Map to physical frame
    Teuchos::Array<double> phys_good_point(spacedim);
    dtk_entity_local_map->mapToPhysicalFrame(dtk_entity, ref_good_point(), phys_good_point());
    for (int d = 0; d < structdim; ++d)
        BOOST_CHECK_CLOSE( good_point[d], phys_good_point[d], percent_tolerance );
    
    // Normal at reference point
    // TODO...
}

