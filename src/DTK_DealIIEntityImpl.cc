#include <no_comment/DTK_DealIIEntityImpl.h>
#include <deal.II/base/geometry_info.h>



// Partial specilization of a function is forbiddent, so we delegue the work to
// some classes that can be specialized
namespace internal
{
  template <int structdim,int dim,int spacedim> struct entity_id;
  template <int dim,int spacedim> struct entity_id<0,dim,spacedim>;

  template <int structdim,int dim,int spacedim>
  struct Entity
  {
    static DataTransferKit::EntityId id(
        Teuchos::RCP<DealIIEntityExtraData<structdim,dim,spacedim>> extra_data)
    {
      auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
            (extra_data)->dealii_tria_accessor;
      DataTransferKit::EntityId entity_id = 0;

      if (structdim == dim) 
      {
        dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
        std::string cell_id = dealii_cell_accessor.id().to_string();
        const unsigned int cell_id_size = cell_id.size();
        for (unsigned int i=0; i<cell_id.size(); ++i)
          entity_id = (static_cast<int>(cell_id[i])<<
              static_cast<int>(std::pow(8,cell_id_size-(i+1)))) | entity_id;
      }
      else
        throw std::runtime_error("entity_id not implemented for faces");

      return entity_id; 
    }
    


    static int ownerRank(Teuchos::RCP<DealIIEntityExtraData<structdim,dim,spacedim>> extra_data)
    {
      auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
        (extra_data)->dealii_tria_accessor;

      if (structdim == dim)
      {
        dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
        return dealii_cell_accessor.subdomain_id();
      }
      else
        throw std::runtime_error("entity_id not implemented for faces");
    }



    static bool inBlock(const int block_id,
        Teuchos::RCP<DealIIEntityExtraData<structdim,dim,spacedim>> extra_data)
    {
      auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
        (extra_data)->dealii_tria_accessor;

      if (structdim == dim)
      {
        dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
        return (dealii_cell_accessor.material_id() == block_id);
      }
      else
        throw std::runtime_error("inBlock not implemented for faces");
    }



    static bool onBoundary(const int boundary_id,
        Teuchos::RCP<DealIIEntityExtraData<structdim,dim,spacedim>> extra_data)
    {
      auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>
        (extra_data)->dealii_tria_accessor;
      if (structdim == dim)
      {
        dealii::CellAccessor<dim,spacedim> dealii_cell_accessor(dealii_tria_accessor);
        // it is on the boundary boundary_id if any of its face is on it
        for (unsigned int face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (dealii_cell_accessor.face(face)->boundary_id() == boundary_id)
            return true;
        }
        return false;
      }
      else
        throw std::runtime_error("onBondary not implemented for faces");
    }
  };

  template <int dim,int spacedim>
  struct Entity<0,dim,spacedim>
  {
    static DataTransferKit::EntityId id(
        Teuchos::RCP<DealIIEntityExtraData<0,dim,spacedim>> extra_data)
    {
      auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<0,dim,spacedim>>
            (extra_data)->dealii_tria_accessor;
      DataTransferKit::EntityId entity_id =
        (*(extra_data->dealii_local_to_global_vertex_id))[dealii_tria_accessor.vertex_index()];

      return entity_id; 
    }



    static int ownerRank(Teuchos::RCP<DealIIEntityExtraData<0,dim,spacedim>> extra_data)
    {
      auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<0,dim,spacedim>>
        (extra_data)->dealii_tria_accessor;
      unsigned int min_rank=-1;
      std::set<typename dealii::Triangulation<dim,spacedim>::active_cell_iterator> adjacent_cells =
        (*(extra_data->dealii_vertex_to_cell))[dealii_tria_accessor.vertex_index()];
      for (auto && adjacent_cell : adjacent_cells)
        if (min_rank > adjacent_cell->subdomain_id())
          min_rank = adjacent_cell->subdomain_id();
      
      return min_rank;
    }



    static bool inBlock(const int block_id,
        Teuchos::RCP<DealIIEntityExtraData<0,dim,spacedim>> extra_data)
    {
      auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<0,dim,spacedim>>
        (extra_data)->dealii_tria_accessor;
      std::set<typename dealii::Triangulation<dim,spacedim>::active_cell_iterator> adjacent_cells =
        (*(extra_data->dealii_vertex_to_cell))[dealii_tria_accessor.vertex_index()];
      for (auto && adjacent_cell : adjacent_cells)
      {
        if (adjacent_cell->material_id() == block_id)
          return true;
      }
      return false;
    }



    static bool onBoundary(const int boundary_id,
        Teuchos::RCP<DealIIEntityExtraData<0,dim,spacedim>> extra_data)
    {
      auto dealii_tria_accessor =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<0,dim,spacedim>>
        (extra_data)->dealii_tria_accessor;
      std::set<typename dealii::Triangulation<dim,spacedim>::active_cell_iterator> adjacent_cells =
        (*(extra_data->dealii_vertex_to_cell))[dealii_tria_accessor.vertex_index()];
      for (auto && adjacent_cell : adjacent_cells)
        for (unsigned int face = 0; face<dealii::GeometryInfo<dim>::faces_per_cell; ++face)
          for (unsigned int v = 0; v<dealii::GeometryInfo<dim>::vertices_per_face; ++v)
            if (adjacent_cell->face(face)->vertex_index(v)==dealii_tria_accessor.vertex_index())
              if (adjacent_cell->face(face)->boundary_id() == boundary_id)
                return true;
      
      return false;
    }
  };
}




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
  return internal::Entity<structdim,dim,spacedim>::id(extra_data);
}



template <int structdim,int dim,int spacedim>
int
DealIIEntityImpl<structdim,dim,spacedim>::
ownerRank() const
{
  return internal::Entity<structdim,dim,spacedim>::ownerRank(extra_data);
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
  return internal::Entity<structdim,dim,spacedim>::inBlock(block_id,extra_data);
}



template <int structdim,int dim,int spacedim>
bool
DealIIEntityImpl<structdim,dim,spacedim>::
onBoundary( const int boundary_id ) const
{ 
  return internal::Entity<structdim,dim,spacedim>::onBoundary(boundary_id,extra_data);
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

