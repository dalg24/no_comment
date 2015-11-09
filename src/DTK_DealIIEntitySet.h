#ifndef DTK_DEALIIENTITYSET_HPP
#define DTK_DEALIIENTITYSET_HPP

#include <no_comment/DTK_DealIIAdjacencies.h>
#include <no_comment/DTK_DealIIEntity.h>
#include <no_comment/DTK_DealIIEntityExtraData.h>
#include <no_comment/DTK_DealIITypes.h>

#include <DTK_EntitySet.hpp>
#include <DTK_Types.hpp>
#include <DTK_Entity.hpp>
#include <DTK_EntityIterator.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit {

template <int dim,int spacedim>
class DealIIEntitySet : public EntitySet
{
public:
    DealIIEntitySet(const Teuchos::RCP<DealIIMesh<dim,spacedim>>& dealii_mesh);

    Teuchos::RCP<const Teuchos::Comm<int>> communicator() const override;

    int physicalDimension() const override;

    void getEntity(const DataTransferKit::EntityId entity_id,
                   const int topological_dimension,
                   DataTransferKit::Entity& entity) const override;

    DataTransferKit::EntityIterator entityIterator(
        const int topoligical_dimension,
        const PredicateFunction& predicate = EntitySet::selectAll) const override;

    void getAdjacentEntities(
        const DataTransferKit::Entity& entity,
        const int adjacent_dimension,
        Teuchos::Array<DataTransferKit::Entity>& adjacent_entities) const override;

    std::string description() const override;

private:
    Teuchos::RCP<DealIIMesh<dim,spacedim>> dealii_triangulation;

    Teuchos::RCP<DealIIAdjacencies<dim,spacedim>> adjacencies;
};

} // end namespace DataTransferKit

#endif
