#include <no_comment/DTK_DealIIEntityImpl.h>
#include <deal.II/base/geometry_info.h>

template <int structdim,int dim,int spacedim>
DealIIEntityImpl<structdim,dim,spacedim>::
DealIIEntityImpl(dealii::TriaAccessor<structdim,dim,spacedim> const & tria_accessor)
: extra_data(Teuchos::rcp(new DealIIEntityExtraData<structdim,dim,spacedim>(tria_accessor)))
{}



template <int structdim,int dim,int spacedim>
DataTransferKit::EntityId
DealIIEntityImpl<structdim,dim,spacedim>::
id() const
{
    auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
            (extra_data)->dealii_tria_accessor;
    DataTransferKit::EntityId entity_id = 0;

    // if entity is a volume element
    if (structdim == dim)
    {
        // Index is unique for only for a given level on a given processor. The
        // level is limited at 11 by p4est.
        unsigned int n_levels = 12;
        entity_id = this->ownerRank() * n_levels * (static_cast<unsigned int>(-1) + 1) +  
            dealii_tria_accessor.level() * (static_cast<unsigned int>(-1) + 1) + 
            dealii_tria_accessor.index();
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
    auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
            (extra_data)->dealii_tria_accessor;
    dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
    return dealii_cell_accessor.subdomain_id();
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
    auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
            (extra_data)->dealii_tria_accessor;
    dealii::Point<spacedim> center = dealii_tria_accessor.center(true);
  
    for (unsigned int i=0; i<spacedim; ++i)
    {
        // Because of round-off the bounding box is made 5% larger than it should be.
        bounds[i] = center[i] - .5*dealii_tria_accessor.extent_in_direction(i);
        bounds[i+3] = center[i] + .5*dealii_tria_accessor.extent_in_direction(i);
    }
}



template <int structdim,int dim,int spacedim>
bool
DealIIEntityImpl<structdim,dim,spacedim>::
inBlock( const int block_id ) const
{        
    auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
            (extra_data)->dealii_tria_accessor;
    if (dim == structdim) {
        dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
        return (dealii_cell_accessor.material_id() == block_id);
    } else {
        throw std::runtime_error("inBlock not implemented for faces and nodes");
    }
}



template <int structdim,int dim,int spacedim>
bool
DealIIEntityImpl<structdim,dim,spacedim>::
onBoundary( const int boundary_id ) const
{ 
    auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
            (extra_data)->dealii_tria_accessor;
    // if the entity is an element
    if (structdim == dim) {
        dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
        // it is on the boundary boundary_id if any of its face is on it
        for (unsigned int face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell; ++face)
            if (dealii_cell_accessor.face(face)->boundary_id() == boundary_id)
                return true;
        return false;
    // if the entity is a face
    } else if (structdim == dim-1) {
        return (dealii_tria_accessor.boundary_id() == boundary_id);
    } else {
        throw std::runtime_error("onBondary not implemented for nodes");
    }
}



template <int structdim,int dim,int spacedim>
Teuchos::RCP<DataTransferKit::EntityExtraData>
DealIIEntityImpl<structdim,dim,spacedim>::
extraData() const
{ 
  return extra_data;
}



template <int structdim,int dim,int spacedim>
void
DealIIEntityImpl<structdim,dim,spacedim>::
describe(
    Teuchos::FancyOStream& os,
    const Teuchos::EVerbosityLevel /*verb_level*/ ) const
{
    os << std::endl;
    os << "---" << std::endl;
    os << "deal.II entity" << std::endl;
    os << "Id: " << id() << std::endl;
    os << "Owner rank: " << ownerRank() << std::endl;
    os << "structdim: " << structdim << std::endl;
    os << "dim: " << dim << std::endl;
    os << "spacedim: " << spacedim << std::endl;
    os << "---" << std::endl;
}
template class DealIIEntityImpl<2,2,2>;
template class DealIIEntityImpl<2,2,3>;
template class DealIIEntityImpl<3,3,3>;

