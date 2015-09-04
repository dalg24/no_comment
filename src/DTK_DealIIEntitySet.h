#ifndef DTK_DEALIIENTITYSET_HPP
#define DTK_DEALIIENTITYSET_HPP

#include "DTK_DealIIAdjacencies.h"
#include "DTK_DealIIEntity.h"
#include "DTK_DealIIEntityExtraData.hpp"

#include <DTK_EntitySet.hpp>
#include <DTK_Types.hpp>
#include <DTK_Entity.hpp>
#include <DTK_EntityIterator.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

template <int dim,int spacedim>
class DealIIEntitySet : public DataTransferKit::EntitySet
{
  public:
    DealIIEntitySet(const Teuchos::RCP<dealii::Triangulation<dim,spacedim>>& dealii_mesh);

    Teuchos::RCP<const Teuchos::Comm<int>> communicator() const override;

    int physicalDimension() const override;

    void getEntity(const DataTransferKit::Entity entity_id,
                   const int topological_dimension,
                   DataTransferKit::Entity& entity) const override;

    DataTransferKit::EntityIterator entityIterator(
        const int topoligical_dimension,
        const std::function<bool(DataTransferKit::Entity)>& predicate) const override;

    void getAdjacentEntities(
        const DataTransferKit::Entity& entity,
        const int adjacent_dimension,
        Teuchos::Array<DataTransferKit::Entity>& adjacent_entities) const override;

    std::string description() const
    { return std::string("deal.II Mesh"); }

  private:
    Teuchos::RCP<dealii::Triangulation> d_dealii_triangulation;

    Teuchos::RCP<DealIIAdjacencies> d_adjacencies;
};

#endif
