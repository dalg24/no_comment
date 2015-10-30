#include "DTK_DealIIEntitySet.h"
#include "DTK_DealIIEntityIterator.h"

#include <Teuchos_DefaultMpiComm.hpp>

template <int dim,int spacedim>
DealIIEntitySet<dim,spacedim>::DealIIEntitySet(
    const Teuchos::RCP<dealii::parallel::distributed::Triangulation<dim,spacedim>> &dealii_mesh)
  :
    d_dealii_triangulation(dealii_mesh),
    d_adjacencies(Teuchos::rcp<DealIIAdjacencies<dim,spacedim>>(
        new DealIIAdjacencies<dim,spacedim>(dealii_mesh)))
{}

template <int dim,int spacedim>
Teuchos::RCP<const Teuchos::Comm<int>>
DealIIEntitySet<dim,spacedim>::
communicator() const
{
    return Teuchos::rcp(
        new Teuchos::MpiComm<int>(d_dealii_triangulation->get_communicator()) );
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
getEntity(const DataTransferKit::EntityId entity_id,
          const int topological_dimension,
          DataTransferKit::Entity& entity) const
{
    if (topological_dimension == dim)
        entity = DealIIEntity<dim,dim,spacedim>(
            *(d_adjacencies->getElemById(entity_id)) );
    else if (topological_dimension == 0)
        entity = DealIIEntity<  0,dim,spacedim>(
            *(d_adjacencies->getNodeById(entity_id)) );
    else
        throw std::runtime_error("not implemented");
}

template <int dim,int spacedim>
DataTransferKit::EntityIterator
DealIIEntitySet<dim,spacedim>::
entityIterator(
    const int topological_dimension,
    const std::function<bool(DataTransferKit::Entity)>& predicate) const
{
    std::ignore = topological_dimension;
    std::ignore = predicate;
    throw std::runtime_error("not implemented");
}


template <int dim,int spacedim>
void
DealIIEntitySet<dim,spacedim>::
getAdjacentEntities(
    const DataTransferKit::Entity& entity,
    const int adjacent_dimension,
    Teuchos::Array<DataTransferKit::Entity>& adjacent_entities) const
{
    if ((entity.physicalDimension() == 0) && (adjacent_dimension  == dim)) {
        std::ignore = adjacent_entities;
        throw std::runtime_error("not implemented");
    } else {
        throw std::runtime_error("not implemented");
    }
}



template class DealIIEntitySet<3,3>;
