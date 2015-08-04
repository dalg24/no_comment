#include <no_comment/DTK_DealIIEntityImpl.h>

template <int structdim,int dim,int spacedim>
DealIIEntityImpl<structdim,dim,spacedim>::
DealIIEntityImpl(std::shared_ptr<dealii::TriaAccessor<structdim,dim,spacedim> const> tria_accessor)
: dealii_tria_accessor(tria_accessor)
{}



template <int structdim,int dim,int spacedim>
DataTransferKit::EntityId
DealIIEntityImpl<structdim,dim,spacedim>::
id() const
{
    DataTransferKit::EntityId entity_id = 0;

    // if entity is a volume element
    if (structdim == dim)
    {
        // Index is unique for only for a given level on a given processor. The
        // level is limited at 11 by p4est.
        unsigned int n_levels = 12;
        entity_id = this->ownerRank() * n_levels * (static_cast<unsigned int>(-1) + 1) +  
          dealii_tria_accessor->level() * (static_cast<unsigned int>(-1) + 1) + 
          dealii_tria_accessor->index();
    }
    else
    {
        throw std::runtime_error("entity_id not implemented for faces and nodes");
    }

    return entity_id; 
}



template <int structdim,int dim,int spacedim>
int
DealIIEntityImpl<structdim,dim,spacedim>::
ownerRank() const
{
    dealii::CellAccessor<dim,spacedim> deallii_cell_accessor(*dealii_tria_accessor);
    return deallii_cell_accessor.subdomain_id();
}



template <int structdim,int dim,int spacedim>
int
DealIIEntityImpl<structdim,dim,spacedim>::
topologicalDimension() const
{
    return dim;
}



template <int structdim,int dim,int spacedim>
int
DealIIEntityImpl<structdim,dim,spacedim>::
physicalDimension() const
{
    return spacedim;
}



template <int structdim,int dim,int spacedim>
void 
DealIIEntityImpl<structdim,dim,spacedim>::
boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    const unsigned int space_dimension = dealii_tria_accessor->space_dimension;
    dealii::Point<space_dimension> center = dealii_tria_accessor->center(true);
  
    for (unsigned int i=0; i<space_dimension; ++i)
    {
        // Because of round-off the bounding box is made 5% larger than it should be.
        bounds[i] = center[i] - .5*dealii_tria_accessor->extent_in_direction(i);
        bounds[i+3] = center[i] + .5*dealii_tria_accessor->extent_in_direction(i);
    }
}



template <int structdim,int dim,int spacedim>
bool
DealIIEntityImpl<structdim,dim,spacedim>::
inBlock( const int block_id ) const
{        
    if (dim == structdim) {
        dealii::CellAccessor<dim,spacedim> deallii_cell_accessor(*dealii_tria_accessor);
        return (deallii_cell_accessor->material_id() == block_id);
    } else {
        throw std::runtime_error("in block not implemented for faces and nodes");
}



//bool
//DealIIEntityImpl::onBoundary( const int boundary_id ) const
//{ 
//  bool on_boundary = false;
//  const dealii::types::boundary_id boundary_indicator = 
//    static_cast<dealii::types::boundary_id>(boundary_id);
//
//  // If tria_accessor is a volume or face, or an edge with dimension equals two,
//  // then we can use the boundary_indicator.
//  if ((tria_accessor->structure_dimension > 1) || 
//      ((tria_accessor->structure_dimension == 1) && (tria_accessor->dimension == 2)))
//    on_boundary = (boundary_indicator == tria_accessor->boundary_indicator()); 
//  else
//  {
//    if (tria_accessor->structure_dimension == 1)
//      for (auto & cell : cell_iterators)
//      {
//        // We need to loop over the faces because the information that we need
//        // does not exist on the edges.
//        for (unsigned int i=0; i<dealii::GeometryInfo<dim>::faces_per_cell; ++i)
//        {
//          for (unsigned int j=0; j<dealii::GeometryInfo<dim>::vertices_per_face; ++j)
//            //TODO need operator == for point and TriaAccessor<0,dim,spacedim>
//            if (cell->face(i)->vertex(j) == tria_accessor)
//              if (cell->face(i)->boundary_indicator == boundary_id)
//                return true;
//        }
//      }
//    else
//    {
//      for (auto & cell : cell_iterators)
//      {
//        // We need to loop over the faces because the information that we need
//        // does not exist on the edges.
//        for (unsigned int i=0; i<dealii::GeometryInfo<dim>::faces_per_cell; ++i)
//        {
//          for (unsigned int j=0; j<dealii::GeometryInfo<dim>::lines_per_face; ++j)
//            if (cell->face(i)->line(j) == tria_accessor)
//              if (cell->face(i)->boundary_indicator == boundary_id)
//                return true;
//        }
//      }
//    }
//  }
//
//  return on_boundary;
//}
//
//
//
//Teuchos::RCP<DataTransferKit::EntityExtraData>
//DealIIEntityImpl::extraData() const
//{ 
//  return Teuchos::rcp(new DataTransferKit::EntityExtraData()); 
//}
