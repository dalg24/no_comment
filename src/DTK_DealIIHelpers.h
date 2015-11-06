#ifndef DTK_DEALIIHELPERS_HPP
#define DTK_DEALIIHELPERS_HPP

#include <no_comment/DTK_DealIIEntityExtraData.h>
#include <DTK_Entity.hpp>
#include <Teuchos_RCP.hpp>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

namespace DataTransferKit {

struct DealIIHelpers
{

    template <int dim,int spacedim>
    static
    dealii::TriaIterator<dealii::CellAccessor<dim,spacedim>>
    getCellIterator(Entity const & entity)
    {
        return dealii::TriaIterator<dealii::CellAccessor<dim,spacedim>>(
            Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<dim,dim,spacedim> const>(
                entity.extraData() )->dealii_tria_accessor );
    }

};

} // end namespace DataTransferKit

#endif
