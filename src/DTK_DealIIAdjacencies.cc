#include <no_comment/DTK_DealIIAdjacencies.h>
#include <no_comment/DTK_DealIIEntity.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/exceptions.h>
#include <boost/mpl/assert.hpp>
#include <type_traits>
#include <limits>
#include <sstream> 

namespace DataTransferKit {

template <int dim,int spacedim>
DealIIAdjacencies<dim,spacedim>::
DealIIAdjacencies(
    Teuchos::RCP<DealIIMesh<dim,spacedim> const> const & tria )
    : tria(tria)
    , local_to_global_vertex_id(dealii::GridTools::compute_local_to_global_vertex_index_map(*tria))
    , vertex_to_cell(dealii::GridTools::vertex_to_cell_map(*tria))
{
    // Reverse the map
    for (auto map_it=local_to_global_vertex_id.cbegin(); 
        map_it!=local_to_global_vertex_id.cend(); ++map_it)
      global_to_local_vertex_id[map_it->second] = map_it->first;


    auto const cell_begin = tria->begin_active();
    auto const cell_end   = tria->end();
    for (auto cell = cell_begin; cell != cell_end; ++cell)
    {
        DealIIElem<dim,spacedim> elem_accessor(tria.get(), cell->level(), cell->index());
        elem_id_to_lvl_index[getId(elem_accessor)] = std::pair<unsigned int, unsigned int>(
            cell->level(), cell->index());
    }
}


template <int dim,int spacedim>
DealIINode<dim,spacedim>
DealIIAdjacencies<dim,spacedim>::
getNodeById(EntityId const id) const
{
    auto ret = global_to_local_vertex_id.find(id);
    AssertThrow(
        ret != global_to_local_vertex_id.end(),
        dealii::ExcMessage("Unknown id") );

    return DealIINode<dim,spacedim>(tria.get(),ret->second);
}


template <int dim,int spacedim>
DealIIElem<dim,spacedim>
DealIIAdjacencies<dim,spacedim>::
getElemById(EntityId const id) const
{
    auto ret = elem_id_to_lvl_index.find(id);
    AssertThrow(
        ret != elem_id_to_lvl_index.end(),
        dealii::ExcMessage("Unknown id") );

    return DealIIElem<dim,spacedim>(tria.get(),ret->second.first,
        ret->second.second);
}


template <int dim,int spacedim>
std::pair<DealIINode<dim,spacedim>,std::vector<dealii::TriaActiveIterator<dealii::CellAccessor<dim,spacedim>>>>
DealIIAdjacencies<dim,spacedim>::
getElemAdjacentToNode(EntityId const id) const
{
    auto ret = global_to_local_vertex_id.find(id);
    AssertThrow(
        ret != global_to_local_vertex_id.end(),
        dealii::ExcMessage("Unknown id") );
  
    DealIINode<dim,spacedim> node(tria.get(),ret->second);
    std::vector<dealii::TriaActiveIterator<dealii::CellAccessor<dim,spacedim>>> adjacent_elem;

    for (auto cell : vertex_to_cell[ret->second])
      adjacent_elem.emplace_back(dealii::TriaIterator<dealii::CellAccessor<dim,spacedim>>(tria.get(),
            cell->level(),cell->index()));


    return std::pair<DealIINode<dim,spacedim>,
      std::vector<dealii::TriaActiveIterator<dealii::CellAccessor<dim,spacedim>>>>(node,adjacent_elem);
}


template <int dim,int spacedim>
EntityId
DealIIAdjacencies<dim,spacedim>::
getId(DealIINode<dim,spacedim> const & node) const
{
    auto ret = local_to_global_vertex_id.find(node.vertex_index());
    if (ret == local_to_global_vertex_id.end())
        throw std::runtime_error("node was not found");
    return ret->second;
}


template <int dim,int spacedim>
EntityId
DealIIAdjacencies<dim,spacedim>::
getId(DealIIElem<dim,spacedim> const & elem) const
{
    EntityId entity_id;
    dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(elem);
    std::string cell_id = dealii_cell_accessor.id().to_string();
    auto colon = cell_id.find(':');
    auto underscore = cell_id.find('_');
    unsigned long long int coarse_cell_id = std::stoul(cell_id.substr(0,underscore));
    unsigned long long int id =  (std::stoul(cell_id.substr(underscore+1,colon-1)) == 0) ? 0 :
      std::stoul(cell_id.substr(colon+1));
    entity_id = coarse_cell_id << std::numeric_limits<unsigned int>::digits | id;

    return entity_id;
}


template <int dim,int spacedim>
DealIIElemIterator<dim,spacedim>
DealIIAdjacencies<dim,spacedim>::
begin_elem() const
{
    return tria->begin_active();
}


template <int dim,int spacedim>
DealIIElemIterator<dim,spacedim>
DealIIAdjacencies<dim,spacedim>::
end_elem() const
{
    return tria->end();
}


template <int dim,int spacedim>
DealIINodeIterator<dim,spacedim>
DealIIAdjacencies<dim,spacedim>::
begin_node() const
{
    return tria->begin_active_vertex();
}


template <int dim,int spacedim>
DealIINodeIterator<dim,spacedim>
DealIIAdjacencies<dim,spacedim>::
end_node() const
{
    return tria->end_vertex();
}


template class DealIIAdjacencies<2,2>;
template class DealIIAdjacencies<2,3>;
template class DealIIAdjacencies<3,3>;

} // end namespace DataTransferKit
