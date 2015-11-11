#include <no_comment/DTK_DealIIEntity.h>
#include <no_comment/DTK_DealIIEntityIterator.h>
#include <type_traits>

namespace DataTransferKit {

// check whether cell is neither locally owned nor ghost
template <int dim,int spacedim>
bool is_artificial(DealIIElemIterator<dim,spacedim> const & tria_iterator)
{
    dealii::CellAccessor<dim,spacedim> cell_accessor(*tria_iterator);
    return ( cell_accessor.subdomain_id() ==
        dealii::numbers::artificial_subdomain_id );
}


template <int structdim,int dim,int spacedim>
DealIIEntityIterator<structdim,dim,spacedim>::
DealIIEntityIterator()
{ }


template <int structdim,int dim,int spacedim>
DealIIEntityIterator<structdim,dim,spacedim>::
DealIIEntityIterator(
    Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies,
    PredicateFunction const &predicate)
    : d_adjacencies(adjacencies)
{
    d_dealii_iterator = adjacencies->begin_elem();
    if (is_artificial(d_dealii_iterator))
        this->operator++();
    d_dealii_iterator_begin = d_dealii_iterator;
    d_dealii_iterator_end = adjacencies->end_elem();
    this->b_predicate = predicate;
}


template <int structdim,int dim,int spacedim>
DealIIEntityIterator<structdim,dim,spacedim>::
DealIIEntityIterator(
    DealIIEntityIterator<structdim,dim,spacedim> const & rhs)
    : EntityIterator()
    , d_dealii_iterator(rhs.d_dealii_iterator)
    , d_dealii_iterator_begin(rhs.d_dealii_iterator_begin)
    , d_dealii_iterator_end(rhs.d_dealii_iterator_end)
    , d_adjacencies(rhs.d_adjacencies)
{
    this->b_predicate = rhs.b_predicate;
}


template <int structdim,int dim,int spacedim>
DealIIEntityIterator<structdim,dim,spacedim>&
DealIIEntityIterator<structdim,dim,spacedim>::
operator=(DealIIEntityIterator<structdim,dim,spacedim> const & rhs)
{
    this->b_predicate = rhs.b_predicate;

    if (&rhs==this)
        return *this;

    d_dealii_iterator = rhs.d_dealii_iterator;
    d_dealii_iterator_begin = rhs.d_dealii_iterator_begin;
    d_dealii_iterator_end = rhs.d_dealii_iterator_end;
    d_adjacencies = rhs.d_adjacencies;

    return *this;
}


template <int structdim,int dim,int spacedim>
EntityIterator &
DealIIEntityIterator<structdim,dim,spacedim>::
operator++()
{
    ++d_dealii_iterator;
    while ((d_dealii_iterator != d_dealii_iterator_end) &&
        is_artificial(d_dealii_iterator))
        ++d_dealii_iterator;

    return *this;
}


template <int structdim,int dim,int spacedim>
Entity&
DealIIEntityIterator<structdim,dim,spacedim>::
operator*(void)
{
    this->operator->();
    return d_current_entity;
}


template <int structdim,int dim,int spacedim>
Entity*
DealIIEntityIterator<structdim,dim,spacedim>::
operator->(void)
{
    d_current_entity = DealIIEntity<structdim,dim,spacedim>(
        *d_dealii_iterator, d_adjacencies );
    return &d_current_entity;
}


template <int structdim,int dim,int spacedim>
bool
DealIIEntityIterator<structdim,dim,spacedim>::
operator==(EntityIterator const &rhs) const
{
    return ( static_cast<DealIIEntityIterator const*>(
        static_cast<DealIIEntityIterator const&>(
            rhs ).b_iterator_impl.get()
        )->d_dealii_iterator == d_dealii_iterator );
}


template <int structdim,int dim,int spacedim>
bool
DealIIEntityIterator<structdim,dim,spacedim>::
operator!=(EntityIterator const & rhs) const
{
  return !(this->operator==(rhs));
}


template <int structdim,int dim,int spacedim>
EntityIterator
DealIIEntityIterator<structdim,dim,spacedim>::
begin() const
{
    DealIIEntityIterator it(
        d_adjacencies,
        this->b_predicate );
    it.d_dealii_iterator = it.d_dealii_iterator_begin;
    return it;
}


template <int structdim,int dim,int spacedim>
EntityIterator
DealIIEntityIterator<structdim,dim,spacedim>::
end() const
{
    DealIIEntityIterator it(
        d_adjacencies,
        this->b_predicate );
    it.d_dealii_iterator = it.d_dealii_iterator_end;
    return it;
}


template <int structdim,int dim,int spacedim>
std::unique_ptr<EntityIterator>
DealIIEntityIterator<structdim,dim,spacedim>::
clone() const
{
    return std::unique_ptr<EntityIterator>(new DealIIEntityIterator(*this));
}


template class DealIIEntityIterator<3,3,3>;

} // end namespace DataTransferKit

