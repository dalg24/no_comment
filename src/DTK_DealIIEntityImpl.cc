#include <no_comment/DTK_DealIIEntityImpl.h>


DealIIEntityImpl::DealIIEntityImpl(TriaAccessor<structdim,dim,spacedim> const &tria_accessor,
    std::vector<CellAccessor<dim,spacedim>> const &cell_accessors) :
  tria_accessor(tria_accessor),
  cell_accessors(cell_accessors)
{}



DataTransferKit::EntityType 
DealIIEntityImpl::entityType() const
{ 
  DataTransferKit::EntityType entity_type;
  switch(tria_accessor->structure_dimension)
  {
    case 0:
      {
        entity_type = DataTransferKit::ENTITY_TYPE_NODE;
        
        break;
      }
    case 1:
      {
        entity_type = DataTransferKit::ENTITY_TYPE_EDGE;
        
        break;
      }
    case 2:
      {
        entity_type = DataTransferKit::ENTITY_TYPE_FACE;
        
        break;
      }
    case 3:
      {
        entity_type = DataTransferKit::ENTITY_TYPE_VOLUME;
        
        break;
      }
    default:
      {
        entity_type = DataTransferKit::ENTITY_TYPE_INVALID;
      }
  }

  return entity_type;
}



DataTransferKit::EntityId
DealIIEntityImpl<dim,spacedim>::id() const
{ 
  // TODO
  DataTransferKit::EntityId entity_id = 0;

  // For cell the ID is unique bur for the other structures, the ID is unique
  // only for a given subdomain. Two different edges on different processors can
  // have the same ID and one edge shared by two different subdomains can have
  // two different IDs.
  if (tria_accessor->structure_dimension == tria_accessor->dimension)
  {
    CellID cell_id = cell_accessors[0]->id();

    // The problem here is that we can't _unpack_ cell_id to compute entity_id.
  }
  else
  {
    entity_id = tria_accessor->index();
  }

  return entity_id; 
}



int
DealIIEntityImpl::ownerRank() const
{ 
  int subddomain_id = -1;

  // TODO 
  // This information only exists for cells
  if (tria_accessor->structure_dimension == tria_accessor->dimension)
    subdomain_id = cell_accessors[0]->subdomain_id(); 

  return subdomain_id;
}



int
DealIIEntityImpl::physicalDimension() const
{ 
  return tria_accessor->space_dimension; 
}



void 
DealIIEntityImpl::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
  const unsigned int space_dimension = tria_accessor->space_dimension;
  Point<space_dimension> center = tria_accessor->barycenter();

  for (unsigned int i=0; i<space_dimension; ++i)
  {
    // To take care of the round-off and the exact definition of the center
    // (should it be the barycenter or should we take into account the
    // manifold?), the bounding is made 10% larger than it should be.
    bounds[i] = center[i] - .55*tria_accessor->extent_in_direction(i);
    bounds[i+3] = center[i] + .55*tria_accessor->extent_in_direction(i);
  }

  for (unsigned int i=space_dimension; i<3; ++i)
  {
    bounds[i] = 0.;
    bounds[i+3] = 0.;
  }
}



bool
DealIIEntityImpl::inBlock( const int block_id ) const
{        
  const types::material_id = static_cast<types::material_id>(block_id);
  // This information only exists for cell. Thus, we loop over the cells that
  // own the tria_accessor.
  for (auto & cell_accessor : cell_accessors)
    if (material_id == cell_accessor->material_id())
      return true

  return false; 
}



bool
DealIIEntityImpl::onBoundary( const int boundary_id ) const
{ 
  return (static_cast<types::boundary_id>(boundary_id) == tria_accessor->boundary_indicator()); 
}



Teuchos::RCP<DataTransferKit::EntityExtraData>
DealIIEntityImpl::extraData() const
{ 
  return Teuchos::rcp(new DataTransferKit::EntityExtraData()); 
}
