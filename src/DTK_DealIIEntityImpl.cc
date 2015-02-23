#include <no_comment/DTK_DealIIEntityImpl.h>

DealIIEntityImpl::
DealIIEntityImpl()
{}

DataTransferKit::EntityType 
DealIIEntityImpl::
entityType() const
{ return DataTransferKit::ENTITY_TYPE_INVALID; }

DataTransferKit::EntityId
DealIIEntityImpl::
id() const
{ return 0; }

int
DealIIEntityImpl::
ownerRank() const
{ return -1; }

int
DealIIEntityImpl::
physicalDimension() const
{ return -1; }

void 
DealIIEntityImpl::
boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{}

bool
DealIIEntityImpl::
inBlock( const int block_id ) const
{ return false; }

bool
DealIIEntityImpl::
onBoundary( const int boundary_id ) const
{ return false; }

Teuchos::RCP<DataTransferKit::EntityExtraData>
DealIIEntityImpl::
extraData() const
{ return Teuchos::rcp(new DataTransferKit::EntityExtraData()); }

