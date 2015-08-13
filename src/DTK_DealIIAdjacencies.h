#ifndef DTK_DEALIIADJACENCIES_HPP
#define DTK_DEALIIADJACENCIES_HPP

#include <no_comment/DTK_DealIIEntity.h>
#include <unordered_map>

template <int dim,int spacedim>
class DealIIAdjacencies
{
public:
    DealIIAdjacencies(std::shared_ptr<dealii::Triangulation<dim,spacedim> const> tria);

    std::shared_ptr<DealIIEntity<dim,dim,spacedim>> getNodeById(DataTransferKit::EntityId const id) const;
    std::shared_ptr<DealIIEntity<0  ,dim,spacedim>> getElemById(DataTransferKit::EntityId const id) const;
private:
    std::unordered_map<
        DataTransferKit::EntityId,
        std::shared_ptr<DealIIEntity<dim,dim,spacedim>>
        > elem_id_map;
    std::unordered_map<
        DataTransferKit::EntityId,
        std::shared_ptr<DealIIEntity<0  ,dim,spacedim>>
        > node_id_map;
    std::unordered_multimap<
        std::shared_ptr<DealIIEntity<0  ,dim,spacedim>>,
        std::shared_ptr<DealIIEntity<dim,dim,spacedim>>
        > node_to_elem_map;
    std::shared_ptr<dealii::Triangulation<dim,spacedim> const> dealii_tria;
};

#endif
