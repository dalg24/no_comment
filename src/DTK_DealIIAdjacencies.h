#ifndef DTK_DEALIIADJACENCIES_HPP
#define DTK_DEALIIADJACENCIES_HPP

#include <no_comment/DTK_DealIIEntity.h>
#include <deal.II/distributed/tria.h>
#include <unordered_map>

namespace DataTransferKit {

template <int dim,int spacedim>
class DealIIAdjacencies
{
using Node = DealIINode<dim,spacedim>*;
using Elem = DealIIElem<dim,spacedim>*;
public:
    DealIIAdjacencies(Teuchos::RCP<DealIIMesh<dim,spacedim> const> dealii_mesh);

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
} // end namespace DataTransferKit

#endif
