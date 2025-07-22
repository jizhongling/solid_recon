cmake_minimum_required(VERSION 3.16)

project(solid_recon_project LANGUAGES CXX)

# Find dependencies
find_package(JANA REQUIRED)
find_package(Threads REQUIRED)
find_package(fmt REQUIRED)

# Compile all sources into executable
file(
  GLOB
  SOURCES
  *.cpp
  *.cc
  *.c
  *.hpp
  *.hh
  *.h)

set(INCLUDE_DIRS ${PROJECT_BINARY_DIR} ${SOLIDRECON_SOURCE_DIR}/src
                 ${PROJECT_SOURCE_DIR} ${JANA_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})
set(LINK_LIBRARIES ${JANA_LIB} ${ROOT_LIBRARIES} ${CMAKE_DL_LIBS}
                   Threads::Threads podio::podio podio::podioRootIO)

# Define executable
add_executable(solid_recon ${SOURCES})
target_include_directories(solid_recon PUBLIC ${INCLUDE_DIRS})
target_link_libraries(solid_recon ${LINK_LIBRARIES})
target_compile_definitions(
  solid_recon PRIVATE SOLIDRECON_APP_VERSION=${CMAKE_PROJECT_VERSION})

# Install executable
install(TARGETS solid_recon DESTINATION bin)
