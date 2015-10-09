#ifndef DTK_DEALIIENTITYITERATOR_H
#define DTK_DEALIIENTITYITERATOR_H
      
#include "no_comment/DTK_DealIIAdjacencies.h"
#include <DTK_EntityIterator.hpp>
#include <DTK_Entity.hpp>
#include <Teuchos_Ptr.hpp>

template <typename DealIIGeomIterator,int dim,int spacedim>
class DealIIEntityIterator : public DataTransferKit::EntityIterator
{
  public:
    DealIIEntityIterator();

    DealIIEntityIterator(
        DealIIGeomIterator dealii_iterator,
        DealIIGeomIterator dealii_iterator_begin,
        DealIIGeomIterator dealii_iterator_end,
        Teuchos::Ptr<dealii::Triangulation<dim,spacedim>> const &dealii_mesh,
        Teuchos::Ptr<DealIIAdjacencies<dim,spacedim>> const &adjacencies,
        PredicateFunction const &predicate);

    DealIIEntityIterator(DealIIEntityIterator<DealIIGeomIterator,dim,spacedim> const &rhs);

    DealIIEntityIterator<DealIIGeomIterator,dim,spacedim> &
      operator= (DealIIEntityIterator<DealIIGeomIterator,dim,spacedim> const &rhs);

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
    Teuchos::Ptr<dealii::Triangulation<dim,spacedim>> d_dealii_mesh;
    Teuchos::Ptr<DealIIAdjacencies<dim,spacedim>> d_adjacencies;
    DataTransferKit::Entity d_current_entity;
};

#include "DTK_DealIIEntityIterator_impl.h"

#endif
