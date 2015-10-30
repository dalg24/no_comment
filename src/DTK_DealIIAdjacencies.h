#ifndef DTK_DEALIIADJACENCIES_HPP
#define DTK_DEALIIADJACENCIES_HPP

#include <no_comment/DTK_DealIIEntity.h>
#include <deal.II/distributed/tria.h>
#include <unordered_map>

template <int dim,int spacedim>
class DealIIAdjacencies
{
using Node = dealii::TriaAccessor<  0,dim,spacedim>*;
using Elem = dealii::TriaAccessor<dim,dim,spacedim>*;
public:
    DealIIAdjacencies(Teuchos::RCP<dealii::parallel::distributed::Triangulation<dim,spacedim> const> tria);

    Node getNodeById(DataTransferKit::EntityId const id) const;
    Elem getElemById(DataTransferKit::EntityId const id) const;
    std::pair<
        typename std::unordered_multimap<Node,Elem>::const_iterator,
        typename std::unordered_multimap<Node,Elem>::const_iterator
        > getElemAdjacentToNode(DataTransferKit::EntityId const id) const;
private:
    std::unordered_map<
        DataTransferKit::EntityId,Elem> elem_id_map;
    std::unordered_map<
        DataTransferKit::EntityId,Node> node_id_map;
    std::unordered_multimap<
        Node,Elem>                      node_to_elem_map;
};

#endif
