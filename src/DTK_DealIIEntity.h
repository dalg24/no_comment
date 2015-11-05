#ifndef DTK_DEALIIENTITY_HPP
#define DTK_DEALIIENTITY_HPP

#include <no_comment/DTK_DealIIAdjacencies.h>
#include <no_comment/DTK_DealIITypes.h>
#include <DTK_Entity.hpp>
#include <Teuchos_RCP.hpp>
#include <map>

namespace DataTransferKit {

template <int structdim,int dim,int spacedim>
class DealIIEntity : public Entity
{
public:
    DealIIEntity(DealIIGeom<structdim,dim,spacedim> const & dealii_tria_accessor,
                 Teuchos::RCP<DealIIAdjacencies<dim,spacedim> const> adjacencies);
};

} // end namespace DataTransferKit

#endif
