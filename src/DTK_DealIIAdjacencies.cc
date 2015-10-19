#include <no_comment/DTK_DealIIAdjacencies.h>
#include <no_comment/DTK_DealIIEntityExtraData.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/base/exceptions.h>
#include <boost/mpl/assert.hpp>
#include <type_traits>
#include <limits>
#include <sstream> 

template <int dim,int spacedim>
DealIIAdjacencies<dim,spacedim>::
DealIIAdjacencies(Teuchos::RCP<dealii::Triangulation<dim,spacedim> const> tria,
    Teuchos::RCP<std::vector<std::set<
    typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>>
    vertex_to_cell,
    Teuchos::RCP<std::map<unsigned int, unsigned long long int>> 
    local_to_global_vertex_id)
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
  auto const cell_begin = tria->begin_active();
  auto const cell_end = tria->end();
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
      node_id_map.emplace((*local_to_global_vertex_id)[vertex_local_id], &(*vertex));
      auto adjacent_cells = (*vertex_to_cell)[vertex_local_id];
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
  AssertThrow(node_id_map.count(id),dealii::ExcMessage("Unknown id"));

  return node_id_map.find(id)->second;
}


template <int dim,int spacedim>
dealii::TriaAccessor<dim,dim,spacedim>*
DealIIAdjacencies<dim,spacedim>::
getElemById(DataTransferKit::EntityId const id) const
{
  AssertThrow(elem_id_map.count(id),dealii::ExcMessage("Unknown id"));

  return elem_id_map.find(id)->second;
}



template class DealIIAdjacencies<2,2>;
template class DealIIAdjacencies<2,3>;
template class DealIIAdjacencies<3,3>;
