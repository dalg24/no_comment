#include <no_comment/DTK_DealIIEntityShapeFunction.h>
#include <no_comment/DTK_DealIIHelpers.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/base/quadrature.h>
#include <DTK_DBC.hpp>

namespace DataTransferKit
{

template <int dim,int spacedim>
DealIIEntityShapeFunction<dim,spacedim>::
DealIIEntityShapeFunction(
    Teuchos::RCP<dealii::DoFHandler<dim,spacedim> const> const & dealii_dof_handler )
    : d_dealii_dof_handler( dealii_dof_handler )
{
//    d_dealii_quadrature =
//        Teuchos::rcp(new dealii::Quadrature<dim>());
//    // @Bruno: should we also somehow pass the mapping?
//    d_dealii_fe_values =
//        Teuchos::rcp(new dealii::FEValues<dim,spacedim>(
//            d_dealii_dof_handler->get_fe(), *d_dealii_quadrature,
//            dealii::update_values | dealii::update_gradients ) );
}



template <int dim,int spacedim>
void
DealIIEntityShapeFunction<dim,spacedim>::
entitySupportIds( 
    const Entity& entity,
    Teuchos::Array<SupportId>& support_ids ) const
{
    auto dof_accessor = *(d_dealii_dof_handler->begin_active());
    dof_accessor.copy_from(
        *DealIIHelpers::getCellIterator<dim,spacedim>(entity) );

    int const dofs_per_cell = d_dealii_dof_handler->get_fe().dofs_per_cell;
    std::vector<dealii::types::global_dof_index> dof_indices(dofs_per_cell);
    dof_accessor.get_dof_indices(dof_indices);

    support_ids.assign(dof_indices.begin(), dof_indices.end());
}



template <int dim,int spacedim>
void
DealIIEntityShapeFunction<dim,spacedim>::
evaluateValue( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<double>& values ) const
{
    dealii::FEValues<dim,spacedim> fe_values(
        d_dealii_dof_handler->get_fe(),
        DealIIHelpers::getPoint<dim>(reference_point),
        dealii::update_values );

    fe_values.reinit(
        DealIIHelpers::getCellIterator<dim,spacedim>(entity) );

    int const dofs_per_cell = fe_values.dofs_per_cell;
    values.resize(dofs_per_cell);
    for (int dof = 0; dof < dofs_per_cell; ++dof)
        values[dof] = fe_values.shape_value(dof, 0);
}



template <int dim,int spacedim>
void
DealIIEntityShapeFunction<dim,spacedim>::
evaluateGradient( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<Teuchos::Array<double> >& gradients ) const
{
    dealii::FEValues<dim,spacedim> fe_values(
        d_dealii_dof_handler->get_fe(),
        DealIIHelpers::getPoint<dim>(reference_point),
        dealii::update_gradients );

    fe_values.reinit(
        DealIIHelpers::getCellIterator<dim,spacedim>(entity) );

    int const dofs_per_cell = fe_values.dofs_per_cell;
    gradients.resize(dofs_per_cell);
    for (int dof = 0; dof < dofs_per_cell; ++dof)
        gradients[dof] = DealIIHelpers::getArray<spacedim>(
            fe_values.shape_grad(dof, 0) );
}

template class DealIIEntityShapeFunction<2,2>;
template class DealIIEntityShapeFunction<3,3>;

} // end namespace DataTransferKit
