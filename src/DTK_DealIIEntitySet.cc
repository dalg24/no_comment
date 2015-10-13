#include "DTK_DealIIEntitySet.h"
#include "DTK_DealIIEntityIterator.h"

#include <Teuchos_DefaultMpiComm.hpp>

template <int dim, int spacedim>
DealIIEntitySet<dim,spacedim>::DealIIEntitySet(
    const Teuchos::RCP<dealii::Triangulation<dim,spacedim>> &dealii_mesh)
  :
    d_dealii_triangulation(dealii_mesh),
    d_adjacencies(Teuchos::rcp<DealIIAdjacencies<dim,spacedim>>(
        new DealIIAdjacencies<dim,spacedim>(dealii_mesh)))
{}
