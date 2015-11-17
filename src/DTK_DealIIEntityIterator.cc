#include <no_comment/DTK_DealIIEntity.h>
#include <no_comment/DTK_DealIIEntityIterator.h>
#include <type_traits>

namespace DataTransferKit {

namespace internal {

template <int structdim,int dim,int spacedim>
struct
EntityIterator
{
    static bool is_artificial(DealIIGeomIterator<structdim,dim,spacedim> const & tria_iterator);
    static DealIIGeomIterator<structdim,dim,spacedim>
    get_begin(Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies);
    static DealIIGeomIterator<structdim,dim,spacedim>
    get_end(Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies);

};

template <int dim,int spacedim>
struct
EntityIterator<dim,dim,spacedim>
{
    static bool is_artificial(
        DealIIGeomIterator<dim,dim,spacedim> const & tria_iterator,
        Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies )
    {
        std::ignore = adjacencies;
        dealii::CellAccessor<dim,spacedim> cell_accessor(*tria_iterator);
        return ( cell_accessor.subdomain_id() ==
            dealii::numbers::artificial_subdomain_id );
    }

    static DealIIGeomIterator<dim,dim,spacedim>
    get_begin(Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies)
    {
        return adjacencies->begin_elem();
    }

    static DealIIGeomIterator<dim,dim,spacedim>
    get_end(Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies)
    {
        return adjacencies->end_elem();
    }

};

template <int dim,int spacedim>
struct
EntityIterator<0,dim,spacedim>
{
    static bool is_artificial(
        DealIIGeomIterator<0,dim,spacedim> const & tria_iterator,
        Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies )
    {
        // TODO: this is a ugly
        // it would be better to have an non throwing function member in
        // adjacencies that can determine if the node exists in the map
        try {
            adjacencies->getId(*tria_iterator);
        } catch (std::runtime_error) {
            return true;
        }
        return false;
    }

    static DealIIGeomIterator<0,dim,spacedim>
    get_begin(Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies)
    {
        return adjacencies->begin_node();
    }

    static DealIIGeomIterator<0,dim,spacedim>
    get_end(Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies)
    {
        return adjacencies->end_node();
    }
};

} // end namespace internal


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
    d_dealii_iterator_end = internal::EntityIterator<structdim,dim,spacedim>::get_end(d_adjacencies);
    d_dealii_iterator = internal::EntityIterator<structdim,dim,spacedim>::get_begin(d_adjacencies);
    if (internal::EntityIterator<structdim,dim,spacedim>::is_artificial(d_dealii_iterator, d_adjacencies))
        this->operator++();
    d_dealii_iterator_begin = d_dealii_iterator;
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
        internal::EntityIterator<structdim,dim,spacedim>::is_artificial(d_dealii_iterator, d_adjacencies))
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


template class DealIIEntityIterator<0,3,3>;
template class DealIIEntityIterator<3,3,3>;

} // end namespace DataTransferKit

