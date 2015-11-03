#include <boost/mpi.hpp>

struct MPIFixture
{
    MPIFixture()
    : argc(boost::unit_test::framework::master_test_suite().argc)
    , argv(boost::unit_test::framework::master_test_suite().argv)
    , env(argc,argv)
    { }
    int argc;
    char **argv;
    boost::mpi::environment env;
    boost::mpi::communicator world;
};
