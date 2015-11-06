#ifndef DTK_DEALIIENTITYEXTRADATA_HPP
#define DTK_DEALIIENTITYEXTRADATA_HPP

#include <no_comment/DTK_DealIITypes.h>
#include <no_comment/DTK_DealIIAdjacencies.h>
#include <DTK_EntityExtraData.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Ptr.hpp>

namespace DataTransferKit {

template <int structdim,int dim,int spacedim>
class DealIIEntityExtraData : public EntityExtraData
{
public:
    DealIIEntityExtraData(DealIIGeom<structdim,dim,spacedim> const & tria_accessor,
        Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> adjacencies)
      : dealii_tria_accessor(tria_accessor)
      , adjacencies(adjacencies)
    { }

    DealIIGeom<structdim,dim,spacedim> const dealii_tria_accessor;
    Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> adjacencies;
};

} // end namespace DataTransferKit

#endif
