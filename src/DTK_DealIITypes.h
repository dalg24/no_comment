#ifndef DTK_DEALIITYPES_HPP
#define DTK_DEALIITYPES_HPP

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

namespace DataTransferKit {

using DealIIVector = dealii::parallel::distributed::Vector<double>;

template <int dim,int spacedim>
using DealIIMesh = dealii::parallel::distributed::Triangulation<dim,spacedim>;

template <int structdim,int dim,int spacedim>
using DealIIGeom = dealii::TriaAccessor<structdim,dim,spacedim>;

template <int structdim,int dim,int spacedim>
using DealIIGeomIterator = dealii::TriaActiveIterator<DealIIGeom<structdim,dim,spacedim>>;

template <int dim,int spacedim>
using DealIINode = DealIIGeom<  0,dim,spacedim>;

template <int dim,int spacedim>
using DealIIElem = DealIIGeom<dim,dim,spacedim>;

template <int dim,int spacedim>
using DealIINodeIterator = DealIIGeomIterator<  0,dim,spacedim>;

template <int dim,int spacedim>
using DealIIElemIterator = DealIIGeomIterator<dim,dim,spacedim>;

}

#endif
