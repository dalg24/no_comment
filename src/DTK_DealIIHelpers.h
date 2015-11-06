#ifndef DTK_DEALIIHELPERS_HPP
#define DTK_DEALIIHELPERS_HPP

#include <no_comment/DTK_DealIIEntityExtraData.h>
#include <DTK_Entity.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/base/point.h>

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

    template <int dim>
    static
    dealii::Point<dim>
    getPoint(Teuchos::ArrayView<double const> const & arr)
    {
        if (arr.size() != static_cast<size_t>(dim))
            std::runtime_error("that's not good");
        dealii::Point<dim> p;
        for (int d = 0; d < dim; ++d)
            p[d] = arr[d];
        return p;
    }

    template <int spacedim>
    static
    Teuchos::Array<double>
    getArray(dealii::Tensor<1,spacedim> const & t)
    {
        Teuchos::Array<double> arr(spacedim);
        for (int d = 0; d < spacedim; ++d)
            arr[d] = t[d];
        return arr;
    }

};

} // end namespace DataTransferKit

#endif
