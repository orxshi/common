cmake_minimum_required(VERSION 3.10.1)

project(common CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})
set(METIS_TEST_RUNS ON)

find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS mpi serialization program_options log thread system filesystem)
find_package(Threads REQUIRED)
find_package(METIS REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wreturn-type")

ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
#ADD_DEFINITIONS(-DNDEBUG)

include_directories(include ${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} ${METIS_INCLUDE_DIRS} ${CGAL_INCLUDE_DIR})
file(GLOB SOURCES "src/*.cpp")
add_library(common SHARED ${SOURCES})
#target_link_libraries(common PUBLIC ${MPI_LIBRARIES} ${Boost_LIBRARIES} ${METIS_LIBRARIES} ${CGAL_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES} Threads::Threads)
target_link_libraries(common PUBLIC ${MPI_LIBRARIES} ${Boost_LIBRARIES} ${METIS_LIBRARIES} Threads::Threads)

if(MPI_COMPILE_FLAGS)
    set_target_properties(common PROPERTIES
        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(common PROPERTIES 
        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
