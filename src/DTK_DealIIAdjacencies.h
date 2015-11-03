#ifndef DTK_DEALIIADJACENCIES_HPP
#define DTK_DEALIIADJACENCIES_HPP

#include <no_comment/DTK_DealIIEntity.h>
#include <deal.II/distributed/tria.h>
#include <unordered_map>

namespace DataTransferKit {

template <int dim,int spacedim>
class DealIIAdjacencies
{
using Node = DealIINode<dim,spacedim>;
using Elem = DealIIElem<dim,spacedim>;
public:
    DealIIAdjacencies(Teuchos::RCP<DealIIMesh<dim,spacedim> const> dealii_mesh);

    Node getNodeById(DataTransferKit::EntityId const id) const;
    Elem getElemById(DataTransferKit::EntityId const id) const;
    std::pair<Node,std::vector<Elem>> getElemAdjacentToNode(DataTransferKit::EntityId const id) const;
private:               
    Teuchos::RCP<dealii::parallel::distributed::Triangulation<dim,spacedim> const> tria;
    std::unordered_map<DataTransferKit::EntityId,unsigned int> global_to_local_vertex_id;
    std::unordered_map<DataTransferKit::EntityId,std::pair<unsigned int,
      unsigned int>> elem_id_to_lvl_index;
    std::vector<std::set<typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>> 
      vertex_to_cell;
};
} // end namespace DataTransferKit

#endif
