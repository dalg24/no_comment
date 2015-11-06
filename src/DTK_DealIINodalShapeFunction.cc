#include <no_comment/DTK_DealIINodalShapeFunction.h>
#include <no_comment/DTK_DealIIHelpers.h>
#include <DTK_DBC.hpp>

namespace DataTransferKit
{

template <int dim,int spacedim>
DealIINodalShapeFunction<dim,spacedim>::
DealIINodalShapeFunction(
    const Teuchos::RCP<dealii::DoFHandler<dim,spacedim>>& dealii_dof_handler )
    : d_dealii_dof_handler( dealii_dof_handler )
{
    d_dealii_quadrature =
        Teuchos::rcp(new dealii::Quadrature<dim>());
    // @Bruno: should we also somehow pass the mapping?
    d_dealii_fe_values =
        Teuchos::rcp(new dealii::FEValues<dim,spacedim>(
            d_dealii_dof_handler->get_fe(), *d_dealii_quadrature,
            dealii::update_values | dealii::update_gradients ) );
}



template <int dim,int spacedim>
void
DealIINodalShapeFunction<dim,spacedim>::
entitySupportIds( 
    const Entity& entity,
    Teuchos::Array<SupportId>& support_ids ) const
{
    auto cell_iterator = DealIIHelpers::getCellIterator<dim,spacedim>(entity);
    
//    // Node case.
//    if ( 0 == entity.topologicalDimension() )
//    {
//	DTK_CHECK( extractGeom<libMesh::Node>(entity)->valid_id() );
//	support_ids.assign( 1, extractGeom<libMesh::Node>(entity)->id() );
//    }
//
//    // Element case.
//    else
//    {
//	Teuchos::Ptr<libMesh::Elem> elem = extractGeom<libMesh::Elem>(entity);
//	int num_nodes = elem->n_nodes();
//	support_ids.resize( num_nodes );
//	for ( int n = 0; n < num_nodes; ++n )
//	{
//	    DTK_CHECK( elem->get_node(n)->valid_id() );
//	    support_ids[n] = elem->get_node(n)->id();
//	}
//    }
}



template <int dim,int spacedim>
void
DealIINodalShapeFunction<dim,spacedim>::
evaluateValue( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<double>& values ) const
{
    auto cell_iterator = DealIIHelpers::getCellIterator<dim,spacedim>(entity);
    d_dealii_fe_values->reinit(cell_iterator);
//    int space_dim = entity.physicalDimension();
//    libMesh::Point lm_reference_point;
//    for ( int d = 0; d < space_dim; ++d )
//    {
//	lm_reference_point(d) = reference_point[d];
//    }
//
//    libMesh::FEComputeData fe_compute_data(
//	d_libmesh_system->get_equation_systems(), lm_reference_point );
//
//    libMesh::FEInterface::compute_data(
//	space_dim,
//	d_libmesh_system->variable_type(0),
//	extractGeom<libMesh::Elem>(entity).getRawPtr(),
//	fe_compute_data );
//
//    values = fe_compute_data.shape;
}



template <int dim,int spacedim>
void
DealIINodalShapeFunction<dim,spacedim>::
evaluateGradient( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<Teuchos::Array<double> >& gradients ) const
{
//    return EntityShapeFunction::evaluateGradient(
//	entity, reference_point, gradients );
}

} // end namespace DataTransferKit
