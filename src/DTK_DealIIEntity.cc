#include <no_comment/DTK_DealIIEntityImpl.h>
#include <no_comment/DTK_DealIIEntity.h>

namespace DataTransferKit {

template <int structdim,int dim,int spacedim>
DealIIEntity<structdim,dim,spacedim>::
DealIIEntity(DealIIGeom<structdim,dim,spacedim> const & dealii_tria_accessor,
             Teuchos::RCP<DealIIAdjacencies<dim,spacedim> const> adjacencies)
{
    this->b_entity_impl = Teuchos::rcp(
        new DealIIEntityImpl<structdim,dim,spacedim>(
            dealii_tria_accessor, adjacencies) );
}

template class DealIIEntity<0,2,2>;
template class DealIIEntity<0,2,3>;
template class DealIIEntity<0,3,3>;
template class DealIIEntity<2,2,2>;
template class DealIIEntity<2,2,3>;
template class DealIIEntity<3,3,3>;

} // end namespace DataTransferKit
