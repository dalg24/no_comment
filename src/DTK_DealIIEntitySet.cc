#include "DTK_DealIIEntitySet.h"

#include <Teuchos_DefaultMpiComm.hpp>

namespace DataTransferKit {

template <int dim,int spacedim>
DealIIEntitySet<dim,spacedim>::
DealIIEntitySet(Teuchos::RCP<DealIIMesh<dim,spacedim>> const & dealii_mesh)
  :
    dealii_triangulation(dealii_mesh),
    adjacencies(Teuchos::rcp<DealIIAdjacencies<dim,spacedim>>(
        new DealIIAdjacencies<dim,spacedim>(dealii_mesh)))
{}

template <int dim,int spacedim>
Teuchos::RCP<const Teuchos::Comm<int>>
DealIIEntitySet<dim,spacedim>::
communicator() const
{
    return Teuchos::rcp(
        new Teuchos::MpiComm<int>(dealii_triangulation->get_communicator()) );
}

template <int dim,int spacedim>
int
DealIIEntitySet<dim,spacedim>::
physicalDimension() const
{
    return spacedim;
}

template <int dim,int spacedim>
void
DealIIEntitySet<dim,spacedim>::
getEntity(const EntityId entity_id,
          const int topological_dimension,
          Entity& entity) const
{
    if (topological_dimension == dim)
        entity = DealIIEntity<dim,dim,spacedim>(
            *(adjacencies->getElemById(entity_id)) );
    else if (topological_dimension == 0)
        entity = DealIIEntity<  0,dim,spacedim>(
            *(adjacencies->getNodeById(entity_id)) );
    else
        throw std::runtime_error("not implemented");
}

template <int dim,int spacedim>
EntityIterator
DealIIEntitySet<dim,spacedim>::
entityIterator(
    const int topological_dimension,
    const std::function<bool(Entity)>& predicate) const
{
    std::ignore = topological_dimension;
    std::ignore = predicate;
    throw std::runtime_error("not implemented");
}


template <int dim,int spacedim>
void
DealIIEntitySet<dim,spacedim>::
getAdjacentEntities(
    const Entity& entity,
    const int adjacent_dimension,
    Teuchos::Array<Entity>& adjacent_entities) const
{
    if ((entity.physicalDimension() == 0) && (adjacent_dimension  == dim)) {
        auto ret = adjacencies->getElemAdjacentToNode(entity.id());
        adjacent_entities.resize(std::distance(ret.first, ret.second));
        auto dtk_entity = adjacent_entities.begin();
        for (auto it = ret.first; it != ret.second; ++it, ++dtk_entity)
            *dtk_entity = DealIIEntity<dim,dim,spacedim>(*(it->second));
        AssertThrow(
            dtk_entity == adjacent_entities.end(),
            dealii::ExcMessage("not good") );
    } else {
        throw std::runtime_error("not implemented");
    }
}

template class DealIIEntitySet<3,3>;

} // end namespace DataTransferKit
