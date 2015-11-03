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
  :
    tria(tria),
    vertex_to_cell(dealii::GridTools::vertex_to_cell_map(*tria))
{
    // checking that dtk has not change its typedef
    static_assert(
        std::is_same<DataTransferKit::EntityId,long unsigned int>::value,
        "dtk entity type is not unsigned long int anymore...");
    // checking that assumptions on integral type sizes are true
    static_assert(
        std::numeric_limits<unsigned long int>::digits ==
        2 * std::numeric_limits<unsigned int>::digits,
        "awwww man...");

    std::map<unsigned int, unsigned long long int> local_to_global_vertex_id =
        dealii::GridTools::compute_local_to_global_vertex_index_map(*tria);
    // Reverse the map
    for (auto map_it=local_to_global_vertex_id.begin(); 
        map_it!=local_to_global_vertex_id.end(); ++map_it)
      global_to_local_vertex_id[map_it->second] = map_it->first;


    auto const cell_begin = tria->begin_active();
    auto const cell_end   = tria->end();
    for (auto cell = cell_begin; cell != cell_end; ++cell)
    {
        dealii::TriaAccessor<dim,dim,spacedim> elem_accessor(tria.get(), cell->level(), cell->index());
        DealIIEntity<dim,dim,spacedim> elem_entity(elem_accessor);
        elem_id_to_lvl_index[elem_entity.id()] = std::pair<unsigned int, unsigned int>(
            cell->level(), cell->index());
    }
}


template <int dim,int spacedim>
dealii::TriaAccessor<0,dim,spacedim>
DealIIAdjacencies<dim,spacedim>::
getNodeById(DataTransferKit::EntityId const id) const
{
    auto ret = global_to_local_vertex_id.find(id);
    AssertThrow(
        ret != global_to_local_vertex_id.end(),
        dealii::ExcMessage("Unknown id") );

    return dealii::TriaAccessor<0,dim,spacedim>(tria.get(),ret->second);
}


template <int dim,int spacedim>
dealii::TriaAccessor<dim,dim,spacedim>
DealIIAdjacencies<dim,spacedim>::
getElemById(DataTransferKit::EntityId const id) const
{
    auto ret = elem_id_to_lvl_index.find(id);
    AssertThrow(
        ret != elem_id_to_lvl_index.end(),
        dealii::ExcMessage("Unknown id") );

    return dealii::TriaAccessor<dim,dim,spacedim>(tria.get(),ret->second.first,
        ret->second.second);
}


template <int dim,int spacedim>
std::pair<dealii::TriaAccessor<0,dim,spacedim>,
  std::vector<dealii::TriaAccessor<dim,dim,spacedim>>>
DealIIAdjacencies<dim,spacedim>::
getElemAdjacentToNode(DataTransferKit::EntityId const id) const
{
    auto ret = global_to_local_vertex_id.find(id);
    AssertThrow(
        ret != global_to_local_vertex_id.end(),
        dealii::ExcMessage("Unknown id") );
  
    dealii::TriaAccessor<0,dim,spacedim> node(tria.get(),ret->second);
    std::vector<dealii::TriaAccessor<dim,dim,spacedim>> adjacent_elem;

    for (auto cell : vertex_to_cell[ret->second])
      adjacent_elem.push_back(dealii::TriaAccessor<dim,dim,spacedim>(tria.get(),
            cell->level(),cell->index()));


    return std::pair<dealii::TriaAccessor<0,dim,spacedim>,
      std::vector<dealii::TriaAccessor<dim,dim,spacedim>>> (node,adjacent_elem);
}


template class DealIIAdjacencies<2,2>;
template class DealIIAdjacencies<2,3>;
template class DealIIAdjacencies<3,3>;
