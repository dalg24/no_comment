FIND_PACKAGE(Trilinos REQUIRED PATHS ${TRILINOS_INSTALL_DIR} NO_DEFAULT_PATH)
INCLUDE_DIRECTORIES(SYSTEM ${DataTransferKit_INCLUDE_DIRS})
MESSAGE("Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("DataTransferKit_INCLUDE_DIRS = ${DataTransferKit_INCLUDE_DIRS}")
