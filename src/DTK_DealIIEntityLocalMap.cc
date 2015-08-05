#include <no_comment/DTK_DealIIEntityLocalMap.h>

template <int dim,int spacedim>
DealIIEntityLocalMap<dim,spacedim>::
DealIIEntityLocalMap(std::shared_ptr<dealii::Mapping<dim,spacedim> const> mapping)
: dealii_mapping(mapping)
{}



template <int dim,int spacedim>
double
DealIIEntityLocalMap<dim,spacedim>::
measure( const DataTransferKit::Entity& entity ) const
{
   return std::nan("");
}



template <int dim,int spacedim>
void
DealIIEntityLocalMap<dim,spacedim>::
centroid(
    const DataTransferKit::Entity& entity,
    const Teuchos::ArrayView<double>& centroid ) const
{
}
