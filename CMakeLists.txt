CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

PROJECT(DealIIDTKAdapters)

#MESSAGE("MPI_INSTALL_DIR=${MPI_INSTALL_DIR}")
#SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${MPI_INSTALL_DIR})

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

INCLUDE(SetupTPLs)

INCLUDE_DIRECTORIES(include)
ADD_LIBRARY(no_comment
    src/DTK_DealIIEntity.cc
    src/DTK_DealIIEntityImpl.cc
    src/DTK_DealIIEntityLocalMap.cc
    src/DTK_DealIIAdjacencies.cc
    src/DTK_DealIIEntityIterator.cc
    src/DTK_DealIIEntitySet.cc
    src/DTK_DealIIEntityShapeFunction.cc
    src/DTK_DealIIEntityIntegrationRule.cc
    src/DTK_DealIIField.cc
    )
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
ENABLE_TESTING()
ADD_SUBDIRECTORY(test)
