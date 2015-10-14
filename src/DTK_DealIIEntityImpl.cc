#include <no_comment/DTK_DealIIEntityImpl.h>
#include <deal.II/base/geometry_info.h>

template <int structdim,int dim,int spacedim>
DealIIEntityImpl<structdim,dim,spacedim>::
DealIIEntityImpl(dealii::TriaAccessor<structdim,dim,spacedim> const & tria_accessor,
    Teuchos::RCP<std::vector<std::set<
    typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>> vertex_to_cell,
    Teuchos::RCP<std::map<unsigned int, unsigned long long int>> local_to_global_vertex_id)
: extra_data(Teuchos::rcp(new DealIIEntityExtraData<structdim,dim,spacedim>(tria_accessor,
        vertex_to_cell, local_to_global_vertex_id)))
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
#if (structdim == dim)
        // Index is unique for only for a given level on a given processor. The
        // level is limited at 11 by p4est.
        unsigned int n_levels = 12;
        entity_id = this->ownerRank() * n_levels * (static_cast<unsigned int>(-1) + 1) +  
            dealii_tria_accessor.level() * (static_cast<unsigned int>(-1) + 1) + 
            dealii_tria_accessor.index();
#elif (structdim == 0)
        entity_id = 
          (*(extra_data->dealii_local_to_global_vertex_id))[dealii_tria_accessor->vertex_index()];
#else
        throw std::runtime_error("entity_id not implemented for faces");
#endif

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
#if (structdim == dim)
    std::cout<<"if (structdim == dim)"<<std::endl;
    std::cout<<"structdim "<<structdim<<std::endl;
    std::cout<<"dim "<<dim<<std::endl;
    std::cout<<"spacedim "<<spacedim<<std::endl;
    dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
    return dealii_cell_accessor.subdomain_id();
#elif (structdim == 0)
    int min_rank=-1;
    std::set<typename dealii::Triangulation<dim,spacedim>::active_cell_iterator> adjacent_cells =
      (*(extra_data->dealii_vertex_to_cell))[dealii_tria_accessor.vertex_index()];
    for (auto && adjacent_cell : adjacent_cells)
    {
      std::cout<<adjacent_cell->center()<<std::endl;
      if (min_rank == -1)
        min_rank = adjacent_cell->subdomain_id();
      else
        if (min_rank > adjacent_cell->subdomain_id())
          min_rank = adjacent_cell->subdomain_id();
    }
    return min_rank;
#else
    throw std::runtime_error("entity_id not implemented for faces");
    return 0;
#endif
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
    // Need to use preprocessor macro otherwise the explicit instantiation
    // creates problem: call vertex_index on a cell

#if (structdim == dim) 
        dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
        return (dealii_cell_accessor.material_id() == block_id);
#elif (structdim == 0)
        std::set<typename dealii::Triangulation<dim,spacedim>::active_cell_iterator> adjacent_cells =
          (*(extra_data->dealii_vertex_to_cell))[dealii_tria_accessor.vertex_index()];
        for (auto && adjacent_cell : adjacent_cells)
        {
          if (adjacent_cell->material_id() == block_id)
            return true;
        }
        return false;
#else 
        throw std::runtime_error("inBlock not implemented for faces");
#endif
}



template <int structdim,int dim,int spacedim>
bool
DealIIEntityImpl<structdim,dim,spacedim>::
onBoundary( const int boundary_id ) const
{ 
    auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
            (extra_data)->dealii_tria_accessor;
    // Need to use preprocessor macro otherwise the explicit instantiation
    // creates problem: call boundary_id on a vertex

    // if the entity is an element
#if (structdim == dim)
        dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
        // it is on the boundary boundary_id if any of its face is on it
        for (unsigned int face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell; ++face)
            if (dealii_cell_accessor.face(face)->boundary_id() == boundary_id)
                return true;
        return false;
    // if the entity is a face
#elif (structdim == dim-1)
        return (dealii_tria_accessor.boundary_id() == boundary_id);
    // if the entity is a vertex
#elif (structdim == 0)
        std::set<typename dealii::Triangulation<dim,spacedim>::active_cell_iterator> adjacent_cells =
          (*(extra_data->dealii_vertex_to_cell))[dealii_tria_accessor.vertex_index()];
        for (auto && adjacent_cell : adjacent_cells)
        {
          for (unsigned int face = 0; face<dealii::GeometryInfo<dim>::faces_per_cell; ++face)
            if (adjacent_cell->face(face)->boundary_id() == boundary_id)
              return true;
        }
        return false;
#else
        throw std::runtime_error("onBondary not implemented for edges");
#endif
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

template class DealIIEntityImpl<0,2,2>;
template class DealIIEntityImpl<0,2,3>;
template class DealIIEntityImpl<0,3,3>;
template class DealIIEntityImpl<2,2,2>;
template class DealIIEntityImpl<2,2,3>;
template class DealIIEntityImpl<3,3,3>;

