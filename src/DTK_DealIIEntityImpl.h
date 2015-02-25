#ifndef DTK_DEALIIENTITYIMPL_HPP
#define DTK_DEALIIENTITYIMPL_HPP

#include <DTK_EntityImpl.hpp>
#include <DTK_EntityExtraData.hpp>
#include <DTK_Types.hpp>

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

using namespace dealii;

// Problem is that tria_accessor for vertices only exist in 1D 
// This probably won't work for 1D mesh.
class DealIIEntityImpl : public DataTransferKit::EntityImpl
{
  public:

    DealIIEntityImpl();

    DataTransferKit::EntityType entityType() const override;

    DataTransferKit::EntityId id() const override;
    
    int ownerRank() const override;

    int physicalDimension() const override;

    void boundingBox( Teuchos::Tuple<double,6>& bounds ) const override;

    bool inBlock( const int block_id ) const override;

    bool onBoundary( const int boundary_id ) const override;

    Teuchos::RCP<DataTransferKit::EntityExtraData> extraData() const override;
};

#endif
