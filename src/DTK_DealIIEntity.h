#ifndef DTK_DEALIIENTITY_HPP
#define DTK_DEALIIENTITY_HPP

#include <DTK_Entity.hpp>
#include <no_comment/DTK_DealIITypes.h>
#include <map>

namespace DataTransferKit {

template <int structdim,int dim,int spacedim>
class DealIIEntity : public Entity
{
public:
    DealIIEntity(DealIIGeom<structdim,dim,spacedim> const & dealii_tria_accessor,
        Teuchos::RCP<std::vector<std::set<
        typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>> 
        vertex_to_cell = Teuchos::ENull::null,
        Teuchos::RCP<std::map<unsigned int, unsigned long long int>> 
        local_to_global_vertex_id = Teuchos::ENull::null);

};

} // end namespace DataTransferKit

#endif
