FIND_PACKAGE(deal.II 8.4 REQUIRED PATHS ${DEAL_II_INSTALL_DIR} NO_DEFAULT_PATH)
INCLUDE(${DEAL_II_TARGET_CONFIG})
MESSAGE("DEAL_II_INCLUDE_DIRS = ${DEAL_II_INCLUDE_DIRS}")
INCLUDE_DIRECTORIES(SYSTEM ${DEAL_II_INCLUDE_DIRS})