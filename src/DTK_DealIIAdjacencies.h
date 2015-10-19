#ifndef DTK_DEALIIADJACENCIES_HPP
#define DTK_DEALIIADJACENCIES_HPP

#include <no_comment/DTK_DealIIEntity.h>
#include <unordered_map>

template <int dim,int spacedim>
class DealIIAdjacencies
{
public:
    DealIIAdjacencies(Teuchos::RCP<dealii::Triangulation<dim,spacedim> const> tria,
    Teuchos::RCP<std::vector<std::set<
    typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>>
    vertex_to_cell,
    Teuchos::RCP<std::map<unsigned int, unsigned long long int>> 
    local_to_global_vertex_id);

    dealii::TriaAccessor<0,dim,spacedim>* getNodeById(DataTransferKit::EntityId const id) const;
    dealii::TriaAccessor<dim,dim,spacedim>* getElemById(DataTransferKit::EntityId const id) const;
private:
    std::unordered_map<
        DataTransferKit::EntityId,
        dealii::TriaAccessor<dim,dim,spacedim>*
        > elem_id_map;
    std::unordered_map<
        DataTransferKit::EntityId,
        dealii::TriaAccessor<0,dim,spacedim>*
        > node_id_map;
    std::unordered_multimap<
        dealii::TriaAccessor<0  ,dim,spacedim>*,
        dealii::TriaAccessor<dim,dim,spacedim>*
        > node_to_elem_map;
};

#endif
