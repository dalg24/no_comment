#### Message Passing Interface (MPI) #########################################
IF(DEFINED MPI_INSTALL_DIR)
    MESSAGE("MPI_INSTALL_DIR=${MPI_INSTALL_DIR}")
    INCLUDE(SetupMPI)
ELSE()
    MESSAGE(FATAL_ERROR "MPI is required. Reconfigure with '-DMPI_INSTALL_DIR=/PATH/TO/MPI'")
ENDIF()

MESSAGE("Setting MPI_CXX_COMPILER as CMAKE_CXX_COMPILER")
SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(--std=c++11 COMPILER_SUPPORTS_CXX11)
IF(NOT COMPILER_SUPPORTS_CXX11)
    MESSAGE(FATAL_ERROR "Compiler has no C++11 support. Please use a different C++ compiler.")
ENDIF()
MESSAGE("Appending -std=c++11 to CMAKE_CXX_FLAGS")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#### deal.II #################################################################
IF(DEFINED DEAL_II_INSTALL_DIR)
    MESSAGE("DEAL_II_INSTALL_DIR=${DEAL_II_INSTALL_DIR}")
    INCLUDE(SetupDealII)
ELSE()
    MESSAGE(FATAL_ERROR "deal.II is required. Reconfigure with '-DDEAL_II_INSTALL_DIR=/PATH/TO/DEAL_II'")
ENDIF()

#### Boost ###################################################################
IF(DEFINED BOOST_INSTALL_DIR)
    MESSAGE("BOOST_INSTALL_DIR=${BOOST_INSTALL_DIR}")
    INCLUDE(SetupBoost)
ELSE()
    MESSAGE(FATAL_ERROR "Boost is required. Reconfigure with '-DBOOST_INSTALL_DIR=/PATH/TO/BOOST'")
ENDIF()

#### Trilinos ################################################################
IF(DEFINED TRILINOS_INSTALL_DIR)
    MESSAGE("TRILINOS_INSTALL_DIR=${TRILINOS_INSTALL_DIR}")
    INCLUDE(SetupTrilinos)
ELSE()
    MESSAGE(FATAL_ERROR "Boost is required. Reconfigure with '-DTRILINOS_INSTALL_DIR=/PATH/TO/TRILINOS'")
ENDIF()

