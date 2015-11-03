#include <no_comment/DTK_DealIIAdjacencies.h>
#include <no_comment/DTK_DealIIEntityExtraData.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/exceptions.h>
#include <boost/mpl/assert.hpp>
#include <type_traits>
#include <limits>
#include <sstream> 

// TODO: move these static assert somewhere else

template <int dim,int spacedim>
DealIIAdjacencies<dim,spacedim>::
DealIIAdjacencies(Teuchos::RCP<dealii::parallel::distributed::Triangulation<dim,spacedim> const> tria)
{
    std::vector<std::set<typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>> vertex_to_cell =
        dealii::GridTools::vertex_to_cell_map(*tria);

    std::map<unsigned int, unsigned long long int> local_to_global_vertex_id =
        dealii::GridTools::compute_local_to_global_vertex_index_map(*tria);

    auto const cell_begin = tria->begin_active();
    auto const cell_end   = tria->end();
    for (auto cell = cell_begin; cell != cell_end; ++cell)
    {
        dealii::TriaAccessor<dim,dim,spacedim> elem_accessor(tria.get(), cell->level(), cell->index());
        DealIIEntity<dim,dim,spacedim> elem_entity(elem_accessor);
        elem_id_map.emplace(elem_entity.id(), &elem_accessor);
        DealIIEntityExtraData<dim,dim,spacedim>* extra_data =
            static_cast<DealIIEntityExtraData<dim,dim,spacedim>*>(elem_entity.extraData().get());
        for (unsigned int i=0; i<dealii::GeometryInfo<dim>::vertices_per_cell; ++i)
        {
            auto vertex = extra_data->dealii_tria_accessor.vertex_iterator(i);
            auto vertex_local_id = extra_data->dealii_tria_accessor.vertex_index(i);
            node_id_map.emplace(local_to_global_vertex_id[vertex_local_id], &(*vertex));
            auto adjacent_cells = vertex_to_cell[vertex_local_id];
            for (auto && adjacent_cell : adjacent_cells)
            {
                dealii::TriaAccessor<dim,dim,spacedim> adjacent_elem(tria.get(),
                    adjacent_cell->level(), adjacent_cell->index());
                node_to_elem_map.emplace(&(*vertex),&adjacent_elem);
            }
        }
    }
}


template <int dim,int spacedim>
dealii::TriaAccessor<0,dim,spacedim>*
DealIIAdjacencies<dim,spacedim>::
getNodeById(DataTransferKit::EntityId const id) const
{
    auto ret = node_id_map.find(id);
    AssertThrow(
        ret != node_id_map.end(),
        dealii::ExcMessage("Unknown id") );

    return ret->second;
}


template <int dim,int spacedim>
dealii::TriaAccessor<dim,dim,spacedim>*
DealIIAdjacencies<dim,spacedim>::
getElemById(DataTransferKit::EntityId const id) const
{
    auto ret = elem_id_map.find(id);
    AssertThrow(
        ret != elem_id_map.end(),
        dealii::ExcMessage("Unknown id") );

    return ret->second;
}


template <int dim,int spacedim>
std::pair<
    typename std::unordered_multimap<
        dealii::TriaAccessor<  0,dim,spacedim>*,
        dealii::TriaAccessor<dim,dim,spacedim>*
        >::const_iterator,
    typename std::unordered_multimap<
        dealii::TriaAccessor<  0,dim,spacedim>*,
        dealii::TriaAccessor<dim,dim,spacedim>*
        >::const_iterator
    >
DealIIAdjacencies<dim,spacedim>::
getElemAdjacentToNode(DataTransferKit::EntityId const id) const
{
      auto ret = node_to_elem_map.equal_range(getNodeById(id));
      AssertThrow(
          ret.first != ret.second,
          dealii::ExcMessage("Unknown id") );

      return ret;
}


template class DealIIAdjacencies<2,2>;
template class DealIIAdjacencies<2,3>;
template class DealIIAdjacencies<3,3>;
