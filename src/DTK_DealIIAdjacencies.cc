#include <no_comment/DTK_DealIIAdjacencies.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/base/exceptions.h>
#include <boost/mpl/assert.hpp>
#include <type_traits>
#include <limits>
#include <sstream> 

template <int dim,int spacedim>
DealIIAdjacencies<dim,spacedim>::
DealIIAdjacencies(std::shared_ptr<dealii::Triangulation<dim,spacedim> const> tria)
: dealii_tria(tria)
{
    // checking that dtk has not change its typedef
    static_assert(
        std::is_same<DataTransferKit::EntityId,long unsigned int>::value,
        "dtk entity type is not unsigned long int anymore...");
    // checking that assumptions on integral type sizes are true
    static_assert(
        std::numeric_limits<unsigned long int>::digits ==
            2 * std::numeric_limits<unsigned int>::digits,
        "awwww man...");
    std::stringstream ss;
    std::string s;
    auto const cell_begin = dealii_tria->begin_active();
    auto const cell_end = dealii_tria->end();
    for (auto cell = cell_begin; cell != cell_end; ++cell)
    {
        auto elem = std::make_shared<DealIIEntity<dim,dim,spacedim>>(*cell);
        ss<<cell->id();
        // convert dealii id into unsigned long int
        std::string dealii_cell_id = ss.str();
        auto colon = dealii_cell_id.find(':');
        auto underscore = dealii_cell_id.find('_');
        unsigned int coarse_cell_id = std::stoul(dealii_cell_id.substr(0, underscore));
        unsigned int id_size = std::stoul(dealii_cell_id.substr(underscore+1,colon-underscore-1));
        long unsigned int id = std::stoul(dealii_cell_id.substr(colon+1));
        std::ignore = id_size;
        if (id > std::numeric_limits<unsigned int>::max())
            throw std::runtime_error("fix me");
        DataTransferKit::EntityId elem_id =  static_cast<long unsigned int>(coarse_cell_id) << 
            std::numeric_limits<unsigned int>::digits | id;

        elem_id_map.emplace(elem_id, elem);

        // empty the string stream
        ss.clear();
        ss.str(std::string());
    }
}


template <int dim,int spacedim>
std::shared_ptr<DealIIEntity<0,dim,spacedim>> 
DealIIAdjacencies<dim,spacedim>::
getNodeById(DataTransferKit::EntityId const id) const
{
  AssertThrow(false,dealii::ExcMessage("Not implemented yet"));

  return node_id_map.find(id)->second;
}


template <int dim,int spacedim>
std::shared_ptr<DealIIEntity<dim,dim,spacedim>> 
DealIIAdjacencies<dim,spacedim>::
getElemById(DataTransferKit::EntityId const id) const
{
  AssertThrow(elem_id_map.count(id),dealii::ExcMessage("Unknown id"));

  return elem_id_map.find(id)->second;
}



template class DealIIAdjacencies<2,2>;
template class DealIIAdjacencies<2,3>;
template class DealIIAdjacencies<3,3>;
