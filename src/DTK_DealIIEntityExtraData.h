#ifndef DTK_DEALIIENTITYEXTRADATA_HPP
#define DTK_DEALIIENTITYEXTRADATA_HPP

#include <DTK_EntityExtraData.hpp>
#include <deal.II/grid/tria_accessor.h>

template <int structdim,int dim,int spacedim>
class DealIIEntityExtraData : public DataTransferKit::EntityExtraData
{
  public:
    DealIIEntityExtraData(dealii::TriaAccessor<structdim,dim,spacedim> const & tria_accessor)
        : dealii_tria_accessor(tria_accessor)
    {}
 
    dealii::TriaAccessor<structdim,dim,spacedim> const dealii_tria_accessor;
};

#endif
