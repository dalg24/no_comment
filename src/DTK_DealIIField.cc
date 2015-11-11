#include <no_comment/DTK_DealIIField.h>

namespace DataTransferKit {

template <int dim,int spacedim>
DealIIField<dim,spacedim>::
DealIIField(
    Teuchos::RCP<dealii::DoFHandler<dim,spacedim> const> const & dof_handler,
    Teuchos::RCP<DealIIVector> const & vector,
    unsigned int const & component )
    : d_dealii_dof_handler(dof_handler)
    , d_dealii_vector(vector)
    , d_component(component)
{
    if (n_components(*d_dealii_dof_handler) != 1)
        throw std::runtime_error("let's restrict us to a single component for now");

    auto const iterator_begin = d_dealii_dof_handler->begin_active();
    auto const iterator_end   = d_dealii_dof_handler->end();
    for (auto it = iterator_begin;
        it != iterator_end; ++it)
    {
        d_support_ids.push_back(-1);
    }
}

template <int dim,int spacedim>
int
DealIIField<dim,spacedim>::
dimension() const
{
    return 1;
//    return d_components.size();
}

template <int dim,int spacedim>
Teuchos::ArrayView<const SupportId>
DealIIField<dim,spacedim>::
getLocalSupportIds() const
{
    return d_support_ids();
}

template <int dim,int spacedim>
double
DealIIField<dim,spacedim>::
readFieldData(
    const SupportId support_id,
    const int dimension ) const
{
    if (dimension != 0)
        throw std::runtime_error("dimension should be zero...");
    dealii::types::global_dof_index const global_dof_index = support_id;
    return (*d_dealii_vector)[global_dof_index];
}


template <int dim,int spacedim>
void
DealIIField<dim,spacedim>::
writeFieldData(
    const SupportId support_id,
    const int dimension,
    const double data)
{
    if (dimension != 0)
        throw std::runtime_error("dimension should be zero...");
    dealii::types::global_dof_index const global_dof_index = support_id;
    (*d_dealii_vector)[global_dof_index] = data;
}


template <int dim,int spacedim>
void
DealIIField<dim,spacedim>::
finalizeAfterWrite()
{
    d_dealii_vector->update_ghost_values();
}

} // end namespace DataTransferKit
