#define BOOST_TEST_MODULE DealIIAdjacencies
#define BOOST_TEST_MAIN
#include <no_comment/DTK_DealIIAdjacencies.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/mpi.h>
#include <boost/test/unit_test.hpp>

struct MPIFixture
{
  MPIFixture() { MPI_Init(NULL,NULL); }
  ~MPIFixture() { MPI_Finalize(); }
};

BOOST_GLOBAL_FIXTURE(MPIFixture);


BOOST_AUTO_TEST_CASE( test_DealIIAdjacencies )
{
  // Probably want to call templated function
  int const dim = 3;
  int const spacedim = 3;

  // Build a mesh
  dealii::parallel::distributed::Triangulation<dim,spacedim> tria_tmp(
      MPI_COMM_WORLD);                                                
  Teuchos::RCP<dealii::parallel::distributed::Triangulation<dim,spacedim>> tria =
    Teuchos::rcpFromRef(tria_tmp);
  dealii::GridGenerator::hyper_rectangle(*tria,
      dealii::Point<spacedim>(-1.0, -2.0, -3.0),
      dealii::Point<spacedim>( 0.0,  0.0,  0.0),
      true);
  tria->refine_global(3);

  std::vector<std::set<
    typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>> vertex_to_cell_tmp(
        dealii::GridTools::vertex_to_cell_map(*tria));
  Teuchos::RCP<std::vector<std::set<
    typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>> vertex_to_cell = 
    Teuchos::rcpFromRef(vertex_to_cell_tmp);
  std::map<unsigned int, unsigned long long int> local_to_global_vertex_id_tmp(
        dealii::GridTools::compute_local_to_global_vertex_index_map(*tria));
  Teuchos::RCP<std::map<unsigned int, unsigned long long int>> local_to_global_vertex_id =
    Teuchos::rcpFromRef(local_to_global_vertex_id_tmp);

  DealIIAdjacencies<dim,spacedim> dealii_adjacencies(tria, vertex_to_cell, 
      local_to_global_vertex_id);
}
