#include <no_comment/DTK_DealIIEntityImpl.h>
#include <no_comment/DTK_DealIIEntity.h>

template <int structdim,int dim,int spacedim>
DealIIEntity<structdim,dim,spacedim>::
DealIIEntity(dealii::TriaAccessor<structdim,dim,spacedim> const & dealii_tria_accessor)
{
    this->b_entity_impl = Teuchos::rcp(
        new DealIIEntityImpl<structdim,dim,spacedim>(
            dealii_tria_accessor) );
}

template class DealIIEntity<3,3,3>;
