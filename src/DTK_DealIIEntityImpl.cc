#include <no_comment/DTK_DealIIEntityImpl.h>


DealIIEntityImpl::DealIIEntityImpl(dealii::TriaAccessor<structdim,dim,spacedim> const &tria_accessor,
    std::vector<active_cell_iterator> const &cell_iterators) :
  tria_accessor(tria_accessor),
  cell_iterators(cell_iterators)
{}



DataTransferKit::EntityType 
DealIIEntityImpl::entityType() const
{ 
  DataTransferKit::EntityType entity_type;
  switch(tria_accessor.structure_dimension)
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
  DataTransferKit::EntityId entity_id = 0;

  if (tria_accessor->structure_dimension == tria_accessor->dimension)
  {
    // Index is unique for only for a given level on a given processor. The
    // level is limited at 11 by p4est.
    unsigned int n_levels = 12;
    entity_id = ownerRank() * n_levels * (static_cast<unsigned int>(-1) + 1) +  
      tria_accessor->level() * (static_cast<unsigned int>(-1) + 1) + 
      tria_accessor->index();
  }
  else
  {
    // Two different edges on different processors can have the same index and 
    // one edge shared by two different subdomains can have two different index.
    entity_id = ownerRank() * (static_cast<unsigned int>(-1) + 1) + tria_accessor->index();
  }

  return entity_id; 
}



int
DealIIEntityImpl::ownerRank() const
{ 
  // This information only exists for cells
  return cell_iterators[0]->subdomain_id(); 
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
  dealii::Point<space_dimension> center = tria_accessor->center(true);

  for (unsigned int i=0; i<space_dimension; ++i)
  {
    // Because of round-off the bounding box is made 5% larger than it should be.
    bounds[i] = center[i] - .525*tria_accessor->extent_in_direction(i);
    bounds[i+3] = center[i] + .525*tria_accessor->extent_in_direction(i);
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
  const dealii::types::material_id = static_cast<dealii::types::material_id>(block_id);
  // This information only exists for cell. Thus, we loop over the cells that
  // own the tria_accessor.
  for (auto & cell : cell_iterators)
    if (material_id == cell->material_id())
      return true

  return false; 
}



bool
DealIIEntityImpl::onBoundary( const int boundary_id ) const
{ 
  bool on_boundary = false;
  const dealii::types::boundary_id boundary_indicator = 
    static_cast<dealii::types::boundary_id>(boundary_id);

  // If tria_accessor is a volume or face, or an edge with dimension equals two,
  // then we can use the boundary_indicator.
  if ((tria_accessor->structure_dimension > 1) || 
      ((tria_accessor->structure_dimension == 1) && (tria_accessor->dimension == 2)))
    on_boundary = (boundary_indicator == tria_accessor->boundary_indicator()); 
  else
  {
    if (tria_accessor->structure_dimension == 1)
      for (auto & cell : cell_iterators)
      {
        // We need to loop over the faces because the information that we need
        // does not exist on the edges.
        for (unsigned int i=0; i<dealii::GeometryInfo<dim>::faces_per_cell; ++i)
        {
          for (unsigned int j=0; j<dealii::GeometryInfo<dim>::vertices_per_face; ++j)
            //TODO need operator == for point and TriaAccessor<0,dim,spacedim>
            if (cell->face(i)->vertex(j) == tria_accessor)
              if (cell->face(i)->boundary_indicator == boundary_id)
                return true;
        }
      }
    else
    {
      for (auto & cell : cell_iterators)
      {
        // We need to loop over the faces because the information that we need
        // does not exist on the edges.
        for (unsigned int i=0; i<dealii::GeometryInfo<dim>::faces_per_cell; ++i)
        {
          for (unsigned int j=0; j<dealii::GeometryInfo<dim>::lines_per_face; ++j)
            if (cell->face(i)->line(j) == tria_accessor)
              if (cell->face(i)->boundary_indicator == boundary_id)
                return true;
        }
      }
    }
  }

  return on_boundary;
}



Teuchos::RCP<DataTransferKit::EntityExtraData>
DealIIEntityImpl::extraData() const
{ 
  return Teuchos::rcp(new DataTransferKit::EntityExtraData()); 
}
