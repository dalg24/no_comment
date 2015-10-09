#ifndef DTK_DEALIIENTITYITERATOR_IMPL_H
#define DTK_DEALIIENTITYITERATOR_IMPL_H

#include <type_traits>
#include <no_comment/DTK_DealIIEntity.h>

template <typename DealIIGeomIterator,int dim,int spacedim>
DealIIEntityIterator<DealIIGeomIterator,dim,spacedim>::DealIIEntityIterator()
{
  this->b_iterator_impl = nullptr;
}


template <typename DealIIGeomIterator,int dim,int spacedim>
DealIIEntityIterator<DealIIGeomIterator,dim,spacedim>::DealIIEntityIterator(
    DealIIGeomIterator dealii_iterator,
    DealIIGeomIterator dealii_iterator_begin,
    DealIIGeomIterator dealii_iterator_end,
    Teuchos::Ptr<dealii::Triangulation<dim,spacedim>> const &dealii_mesh,
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


template <typename DealIIGeomIterator,int dim,int spacedim>
DealIIEntityIterator<DealIIGeomIterator,dim,spacedim>&
DealIIEntityIterator<DealIIGeomIterator,dim,spacedim>::operator=(
    DealIIEntityIterator<DealIIGeomIterator,dim,spacedim> const & rhs)
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


template <typename DealIIGeomIterator,int dim,int spacedim>
DataTransferKit::EntityIterator &
DealIIEntityIterator<DealIIGeomIterator,dim,spacedim>::operator++()
{
  ++d_dealii_iterator;

  return *this;
}


template <typename DealIIGeomIterator,int dim,int spacedim>
DataTransferKit::Entity*
DealIIEntityIterator<DealIIGeomIterator,dim,spacedim>::operator->(void)
{
  d_current_entity = DealIIEntity<
    typename std::remove_pointer<
    typename DealIIGeomIterator::value_type> (
        Teuchos::ptr(*d_dealii_iterator), d_dealii_mesh, d_adjacencies);

  return &d_current_entity;
}


template <typename DealIIGeomIterator,int dim,int spacedim>
bool DealIIEntityIterator<DealIIGeomIterator,dim,spacedim>::operator==(
    DataTransferKit::EntityIterator const &)
{
  DealIIEntityIterator const* rhs_it = 
    static_cast<DealIIEntityIterator const*>(&rhs);
  DealIIEntityIterator const* rhs_it_impl = 
    static_cast<DealIIEntityIterator const*>(rhs_it->b_iterator_impl);

  return (rhs_it_impl->d_dealii_iterator == d_dealii_iterator);
}


template <typename DealIIGeomIterator,int dim,int spacedim>
bool DealIIEntityIterator<DealIIGeomIterator,dim,spacedim>::operator!=(
    DataTransferKit::EntityIterator & const rhs) const
{
  DealIIEntityIterator const* rhs_it =
    static_cast<DealIIEntityIterator const*>(&rhs);
  DealIIEntityIterator const* rhs_it_impl =
    static_cast<DealIIEntityIterator const*>(rhs_it->b_iterator_impl);

  return (rhs_it_impl->d_dealii_iterator != d_dealii_iterator);
}


template <typename DealIIGeomIterator,int dim,int spacedim>
DataTransferKit::EntityIterator<DealIIGeomIterator,dim,spacedim>::begin() const
{
  return DealIIEntityIterator(d_dealii_iterator_begin,
      d_dealii_iterator_begin,
      d_dealii_iterator_end,
      d_dealii_mesh,
      d_adjacencies,
      this->b_predicate);
}



template <typename DealIIGeomIterator,int dim,int spacedim>
DataTransferKit::EntityIterator
DealIIEntityIterator<DealIIGeomIterator,dim,spacedim>::end() const
{
  return DealIIEntityIterator(d_dealii_iterator_end,
      d_dealii_iterator_begin,
      d_dealii_iterator_end,
      d_dealii_mesh,
      d_adjacencies,
      this->b_predicate);
}



template <typename DealIIGeomIterator,int dim,int spacedim>
DataTransferKit::EntityIterator*
DealIIEntityIterator<DealIIGeomIterator,int dim,int spacedim>::clone() const
{
  return new DealIIEntityIterator(*this);
}

#endif
