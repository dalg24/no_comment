#ifndef DTK_DEALIIENTITYITERATOR_H
#define DTK_DEALIIENTITYITERATOR_H
      
#include <no_comment/DTK_DealIIAdjacencies.h>
#include <DTK_EntityIterator.hpp>
#include <DTK_Entity.hpp>
#include <Teuchos_Ptr.hpp>

namespace DataTransferKit {

template <int structdim,int dim,int spacedim>
class DealIIEntityIterator : public EntityIterator
{
  public:
    DealIIEntityIterator();

    DealIIEntityIterator(
        Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> const & adjacencies,
        PredicateFunction const & predicate );

    DealIIEntityIterator(DealIIEntityIterator<structdim,dim,spacedim> const &rhs);

    DealIIEntityIterator<structdim,dim,spacedim> &
    operator=(DealIIEntityIterator<structdim,dim,spacedim> const &rhs);

    EntityIterator& operator++() override;

    Entity& operator*(void) override;

    Entity* operator->(void) override;

    bool operator==(EntityIterator const &rhs) const override;

    bool operator!=(EntityIterator const &rhs) const override;

    EntityIterator begin() const override;

    EntityIterator end() const override;

    std::unique_ptr<EntityIterator> clone() const override;

  private:
    DealIIGeomIterator<structdim,dim,spacedim> d_dealii_iterator;
    DealIIGeomIterator<structdim,dim,spacedim> d_dealii_iterator_begin;
    DealIIGeomIterator<structdim,dim,spacedim> d_dealii_iterator_end;
    Teuchos::Ptr<DealIIAdjacencies<dim,spacedim> const> d_adjacencies;
    Entity d_current_entity;
};

} // end namespace DataTransferKit

#endif
