struct MPIFixture
{
    MPIFixture()
    : argc(boost::unit_test::framework::master_test_suite().argc)
    , argv(boost::unit_test::framework::master_test_suite().argv)
    {
        MPI_Init(&argc,&argv);
    }

    ~MPIFixture()
    {
        MPI_Finalize();
    }

    int argc;
    char **argv;
};
