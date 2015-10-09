#include "DTK_DealIIEntitySet.h"
#include "DTK_DealIIEntityIterator.h"

#include <Teuchos_DefaultMpiComm.hpp>

template <int dim, int spacedim>
DTK_DealIIEntitySet<dim,spacedim>::DTK_DealIIEntitySet(
    const Teuchos::RCP<dealii::Triangulation<dim,spacedim> &dealii_mesh)
  :
    d_dealii_triangulation(dealii_mesh),
    d_adjancencies(new DealIIAdjacencies(dealii_mesh))
{}
