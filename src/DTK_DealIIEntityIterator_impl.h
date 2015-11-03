#ifndef DTK_DEALIIENTITYITERATOR_IMPL_H
#define DTK_DEALIIENTITYITERATOR_IMPL_H

#include <type_traits>
#include <no_comment/DTK_DealIIEntity.h>

namespace DataTransferKit {

template <int structdim,int dim,int spacedim>
DealIIEntityIterator<structdim,dim,spacedim>::
DealIIEntityIterator()
{
  this->b_iterator_impl = nullptr;
}


template <int structdim,int dim,int spacedim>
DealIIEntityIterator<structdim,dim,spacedim>::
DealIIEntityIterator(
    DealIIGeomIterator<structdim,dim,spacedim> dealii_iterator,
    DealIIGeomIterator<structdim,dim,spacedim> dealii_iterator_begin,
    DealIIGeomIterator<structdim,dim,spacedim> dealii_iterator_end,
    Teuchos::Ptr<dealii::parallel::distributed::Triangulation<dim,spacedim>> const &dealii_mesh,
    Teuchos::Ptr<DealIIAdjacencies<dim,spacedim>> const &adjacencies,
    PredicateFunction const &predicate)
  :
    d_dealii_iterator(dealii_iterator),
    d_dealii_iterator_begin(dealii_iterator_begin),
    d_dealii_iterator_end(dealii_iterator_end),
    d_dealii_mesh(dealii_mesh),
    d_adjacencies(adjacencies)
{
  this->b_iterator_impl = nullptr;
  this->b_predicate = predicate;
}


template <int structdim,int dim,int spacedim>
DealIIEntityIterator<structdim,dim,spacedim>&
DealIIEntityIterator<structdim,dim,spacedim>::
operator=(DealIIEntityIterator<structdim,dim,spacedim> const & rhs)
{
  this->b_iterator_impl = nullptr;
  this->b_predicate = rhs.b_predicate;
  if (&rhs==this)
    return *this;
  d_dealii_iterator = rhs.d_dealii_iterator;
  d_dealii_iterator_begin = rhs.d_dealii_iterator_begin;
  d_dealii_iterator_end = rhs.d_dealii_iterator_end;
  d_dealii_mesh = rhs.d_dealii_mesh;
  d_adjacencies = rhs.d_adjacencies;

  return *this;
}


template <int structdim,int dim,int spacedim>
EntityIterator &
DealIIEntityIterator<structdim,dim,spacedim>::
operator++()
{
  ++d_dealii_iterator;

  return *this;
}


template <int structdim,int dim,int spacedim>
Entity*
DealIIEntityIterator<structdim,dim,spacedim>::
operator->(void)
{
  // This will probably not work but it compiles which is the most important :-)
  d_current_entity = DealIIEntity<
         structdim,dim,spacedim
    >( Teuchos::ptr(*d_dealii_iterator), d_dealii_mesh, d_adjacencies);

  return &d_current_entity;
}


template <int structdim,int dim,int spacedim>
bool DealIIEntityIterator<structdim,dim,spacedim>::
operator==(EntityIterator const &rhs) const
{
  DealIIEntityIterator const* rhs_it = 
    static_cast<DealIIEntityIterator const*>(&rhs);
  DealIIEntityIterator const* rhs_it_impl = 
    static_cast<DealIIEntityIterator const*>(rhs_it->b_iterator_impl);

  return (rhs_it_impl->d_dealii_iterator == d_dealii_iterator);
}


template <int structdim,int dim,int spacedim>
bool DealIIEntityIterator<structdim,dim,spacedim>::
operator!=(EntityIterator const & rhs) const
{
  DealIIEntityIterator const* rhs_it =
    static_cast<DealIIEntityIterator const*>(&rhs);
  DealIIEntityIterator const* rhs_it_impl =
    static_cast<DealIIEntityIterator const*>(rhs_it->b_iterator_impl);

  return (rhs_it_impl->d_dealii_iterator != d_dealii_iterator);
}


template <int structdim,int dim,int spacedim>
EntityIterator
DealIIEntityIterator<structdim,dim,spacedim>::
begin() const
{
  return DealIIEntityIterator(d_dealii_iterator_begin,
      d_dealii_iterator_begin,
      d_dealii_iterator_end,
      d_dealii_mesh,
      d_adjacencies,
      this->b_predicate);
}



template <int structdim,int dim,int spacedim>
EntityIterator
DealIIEntityIterator<structdim,dim,spacedim>::
end() const
{
  return DealIIEntityIterator(d_dealii_iterator_end,
      d_dealii_iterator_begin,
      d_dealii_iterator_end,
      d_dealii_mesh,
      d_adjacencies,
      this->b_predicate);
}



template <int structdim,int dim,int spacedim>
EntityIterator*
DealIIEntityIterator<structdim,dim,spacedim>::
clone() const
{
  return new DealIIEntityIterator(*this);
}

} // end namespace DataTransferKit

#endif
