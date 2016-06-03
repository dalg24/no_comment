#ifndef DTK_DEALIIADJACENCIES_HPP
#define DTK_DEALIIADJACENCIES_HPP

#include <no_comment/DTK_DealIITypes.h>
#include <deal.II/distributed/tria.h>
#include <DTK_Types.hpp>
#include <Teuchos_RCP.hpp>
#include <unordered_map>

namespace DataTransferKit {

template <int dim,int spacedim>
class DealIIAdjacencies
{
using Node = DealIINode<dim,spacedim>;
using Elem = DealIIElem<dim,spacedim>;
using CellIterator = dealii::TriaActiveIterator<dealii::CellAccessor<dim,spacedim>>;
using Mesh = DealIIMesh<dim,spacedim>;
public:
    DealIIAdjacencies(Teuchos::RCP<Mesh const> const & dealii_mesh);

    Node getNodeById(EntityId const id) const;
    Elem getElemById(EntityId const id) const;
//    std::pair<Node,std::vector<Elem>> getElemAdjacentToNode(EntityId const id) const;
    std::pair<Node,std::vector<CellIterator>> getElemAdjacentToNode(EntityId const id) const;
    EntityId getId(Node const & node) const;
    EntityId getId(Elem const & elem) const;
    DealIIElemIterator<dim,spacedim> begin_elem() const;
    DealIIElemIterator<dim,spacedim> end_elem() const;
    DealIINodeIterator<dim,spacedim> begin_node() const;
    DealIINodeIterator<dim,spacedim> end_node() const;
private:               
    Teuchos::RCP<Mesh const> tria;
    std::map<unsigned int,unsigned long long int> local_to_global_vertex_id;
    std::unordered_map<EntityId,unsigned int> global_to_local_vertex_id;
    std::unordered_map<EntityId,std::pair<unsigned int,
      unsigned int>> elem_id_to_lvl_index;
    std::vector<std::set<CellIterator>>
      vertex_to_cell;
};
} // end namespace DataTransferKit

#endif
