#ifndef DTK_DEALIIENTITYEXTRADATA_HPP
#define DTK_DEALIIENTITYEXTRADATA_HPP

#include <DTK_EntityExtraData.hpp>
#include <Teuchos_ArrayView.hpp>
#include <deal.II/grid/tria_accessor.h>

template <int structdim,int dim,int spacedim>
class DealIIEntityExtraData : public DataTransferKit::EntityExtraData
{
  public:
    DealIIEntityExtraData(dealii::TriaAccessor<structdim,dim,spacedim> const & tria_accessor,
        Teuchos::RCP<std::vector<std::set<
        typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>> vertex_to_cell,
        Teuchos::RCP<std::map<unsigned int, unsigned long long int>> local_to_global_vertex_id)
      : 
        dealii_tria_accessor(tria_accessor),
        dealii_vertex_to_cell(vertex_to_cell),
        dealii_local_to_global_vertex_id(local_to_global_vertex_id)
  {}

    dealii::TriaAccessor<structdim,dim,spacedim> const dealii_tria_accessor;
    // Two pointers per vertex is a lot of pointers :/ It should only be necessary
    // to store one pointer per processor.
    Teuchos::RCP<std::vector<std::set<
      typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>> dealii_vertex_to_cell;
    Teuchos::RCP<std::map<unsigned int, unsigned long long int>> dealii_local_to_global_vertex_id;
};

#endif
