#ifndef DTK_DEALIIENTITY_HPP
#define DTK_DEALIIENTITY_HPP

#include <DTK_Entity.hpp>
#include <deal.II/grid/tria_accessor.h>

template <int structdim,int dim,int spacedim>
class DealIIEntity : public DataTransferKit::Entity
{
public:
    DealIIEntity(dealii::TriaAccessor<structdim,dim,spacedim> const & dealii_tria_accessor);

};

#endif
