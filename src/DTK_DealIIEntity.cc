#include <no_comment/DTK_DealIIEntityImpl.h>
#include <no_comment/DTK_DealIIEntity.h>

namespace DataTransferKit {

template <int structdim,int dim,int spacedim>
DealIIEntity<structdim,dim,spacedim>::
DealIIEntity(DealIIGeom<structdim,dim,spacedim> const & dealii_tria_accessor,
    Teuchos::RCP<std::vector<std::set<
    typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>> vertex_to_cell,
    Teuchos::RCP<std::map<unsigned int, unsigned long long int>> local_to_global_vertex_id)
{
    this->b_entity_impl = Teuchos::rcp(
        new DealIIEntityImpl<structdim,dim,spacedim>(
            dealii_tria_accessor, vertex_to_cell, local_to_global_vertex_id) );
}

template class DealIIEntity<0,2,2>;
template class DealIIEntity<0,2,3>;
template class DealIIEntity<0,3,3>;
template class DealIIEntity<2,2,2>;
template class DealIIEntity<2,2,3>;
template class DealIIEntity<3,3,3>;

} // end namespace DataTransferKit
