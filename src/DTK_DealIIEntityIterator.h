#ifndef DTK_DEALIIENTITYITERATOR_H
#define DTK_DEALIIENTITYITERATOR_H
      
#include "no_comment/DTK_DealIIAdjacencies.h"
#include <DTK_EntityIterator.hpp>
#include <DTK_Entity.hpp>
#include <Teuchos_Ptr.hpp>

template <int structdim,int dim,int spacedim>
class DealIIEntityIterator : public DataTransferKit::EntityIterator
{
  using DealIIGeomIterator =
      dealii::TriaActiveIterator<dealii::TriaAccessor<structdim,dim,spacedim>>;
  public:
    DealIIEntityIterator();

    DealIIEntityIterator(
        DealIIGeomIterator dealii_iterator,
        DealIIGeomIterator dealii_iterator_begin,
        DealIIGeomIterator dealii_iterator_end,
        Teuchos::Ptr<dealii::parallel::distributed::Triangulation<dim,spacedim>> const &dealii_mesh,
        Teuchos::Ptr<DealIIAdjacencies<dim,spacedim>> const &adjacencies,
        DataTransferKit::PredicateFunction const &predicate);

    DealIIEntityIterator(DealIIEntityIterator<structdim,dim,spacedim> const &rhs);

    DealIIEntityIterator<structdim,dim,spacedim> &
      operator= (DealIIEntityIterator<structdim,dim,spacedim> const &rhs);

    ~DealIIEntityIterator();

    DataTransferKit::EntityIterator& operator++() override;

    DataTransferKit::Entity& operator*(void) override;

    DataTransferKit::Entity* operator->(void) override;

    bool operator== (DataTransferKit::EntityIterator const &rhs) const override;

    bool operator!= (DataTransferKit::EntityIterator const &rhs) const override;

    DataTransferKit::EntityIterator begin() const override;

    DataTransferKit::EntityIterator end() const override;

    DataTransferKit::EntityIterator* clone() const override;

  private:
    DealIIGeomIterator d_dealii_iterator;
    DealIIGeomIterator d_dealii_iterator_begin;
    DealIIGeomIterator d_dealii_iterator_end;
    Teuchos::Ptr<dealii::parallel::distributed::Triangulation<dim,spacedim>> d_dealii_mesh;
    Teuchos::Ptr<DealIIAdjacencies<dim,spacedim>> d_adjacencies;
    DataTransferKit::Entity d_current_entity;
};

#include "DTK_DealIIEntityIterator_impl.h"

#endif
