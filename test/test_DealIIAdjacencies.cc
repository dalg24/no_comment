#define BOOST_TEST_MODULE DealIIAdjacencies
#define BOOST_TEST_MAIN
#include <no_comment/DTK_DealIIAdjacencies.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/mpi.h>
#include <boost/test/unit_test.hpp>

#include "MPIFixture.cc"

BOOST_FIXTURE_TEST_CASE( test_DealIIAdjacencies, MPIFixture )
{
  // Probably want to call templated function
  int const dim = 3;
  int const spacedim = 3;

  // Build a mesh
  Teuchos::RCP<DataTransferKit::DealIIMesh<dim,spacedim>> tria =
    Teuchos::rcp(new DataTransferKit::DealIIMesh<dim,spacedim>(world));

  dealii::GridGenerator::hyper_rectangle(*tria,
      dealii::Point<spacedim>(-1.0, -2.0, -3.0),
      dealii::Point<spacedim>( 0.0,  0.0,  0.0),
      true);
  tria->refine_global(3);

  DataTransferKit::DealIIAdjacencies<dim,spacedim> dealii_adjacencies(tria);

  double const tolerance = 1.0e-15;
  dealii::TriaAccessor<0,dim,spacedim> node = dealii_adjacencies.getNodeById(0);
  BOOST_CHECK_SMALL(node.vertex()[0]+1.,tolerance);
  BOOST_CHECK_SMALL(node.vertex()[1]+2.,tolerance);
  BOOST_CHECK_SMALL(node.vertex()[2]+3.,tolerance);
  node = dealii_adjacencies.getNodeById(705);
  BOOST_CHECK_SMALL(node.vertex()[0]+0.125,tolerance);
  BOOST_CHECK_SMALL(node.vertex()[1]+0.75,tolerance);
  BOOST_CHECK_SMALL(node.vertex()[2]+0.375,tolerance);

  unsigned int level = 3;
  unsigned int index = 2;
  dealii::TriaAccessor<dim,dim,spacedim> elem_accessor(tria.get(), level, index);
  DataTransferKit::DealIIEntity<dim,dim,spacedim> elem_entity(elem_accessor);
  dealii::TriaAccessor<dim,dim,spacedim> elem = dealii_adjacencies.getElemById(
    elem_entity.id());
  BOOST_CHECK_EQUAL(elem.level(),level);
  BOOST_CHECK_EQUAL(elem.index(),index);


  std::pair<DataTransferKit::DealIINode<dim,spacedim>, std::vector<
    DataTransferKit::DealIIElem<dim,spacedim>>> node_elems =
    dealii_adjacencies.getElemAdjacentToNode(705);
  BOOST_CHECK_SMALL(node_elems.first.vertex()[0]+0.125,tolerance);
  BOOST_CHECK_SMALL(node_elems.first.vertex()[1]+0.75,tolerance);
  BOOST_CHECK_SMALL(node_elems.first.vertex()[2]+0.375,tolerance);

  index = 488;
  for (auto adjacent_elem : node_elems.second)
  {
    BOOST_CHECK_EQUAL(adjacent_elem.level(),level);
    BOOST_CHECK_EQUAL(adjacent_elem.index(),index);
    ++index;
  }
}
