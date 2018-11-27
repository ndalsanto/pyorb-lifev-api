CMAKE_MINIMUM_REQUIRED(VERSION 3.9)

project(lifev-pyorb VERSION 1.0.0 DESCRIPTION "lifev-pyorb description")

# required for unit testing using CMake's ctest command
include(CTest)
enable_testing()

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(ExternalProject)

# set(CMAKE_BUILD_TYPE Release)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

SET( INCLUDED_HEADERS
    ${MPI_INCLUDE_DIRS}
    ${Trilinos_INCLUDE_DIRS}
    ${TPL_LifeV_INCLUDE_DIRS} )

SET( LINKED_LIBRARIES_DIRS
   ${MPI_LIBRARY_DIRS}
   ${MPI_LIBRARY_DIRS}
   ${Trilinos_LIBRARY_DIRS}
   ${TPL_LifeV_LIBRARY_DIRS}
    )

SET( LINKED_LIBRARIES
     libmpi.so
     libmpi_cxx.so
     liblifevcore.so
     liblifeveta.so
      )

include_directories( ${INCLUDED_HEADERS} )
link_directories( ${LINKED_LIBRARIES_DIRS} )

add_subdirectory(include)
add_subdirectory(sources)

set(HEADER
    ${HEADER_FILES}
  )

set(SOURCE
    ${SOURCE_FILES}
  )

add_library(tessst SHARED ${SOURCE} )
target_link_libraries(tessst ${LINKED_LIBRARIES} )


add_subdirectory(test)


# SET(LIFEV_LIBRARIES_PATHS ${TPL_LifeV_LIBRARY_DIRS}/lib/liblifevcore.so )
