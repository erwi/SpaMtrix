CMAKE_MINIMUM_REQUIRED (VERSION 3.0.0)
PROJECT (SpaMtrix)
# let user define whether to compile with support for OpenMP
# value should be ON/OFF
# the definition USES_OPENMP is added, if enabled and can be
# checked for in code using #ifdef USES_OPENMP
option (USES_OPENMP "whether to use OpenMP?" OFF)
if (USES_OPENMP)
    # CHECK WHETHER OPENMP IS AVAILABLE
    find_package(OpenMP)
    if (OPENMP_FOUND)
        # TODO: get rid of these variables
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${OpenMP_CXX_FLAGS}") #-lpthread -lgomp")
        add_definitions(-DUSES_OPENMP)
    endif(OPENMP_FOUND)
endif(USES_OPENMP)
#
# add build directories
add_subdirectory(lib)
add_subdirectory(examples)

# so that tests are only run if this is the parent level project, not
# of SpaMtrix is included as a child project in another project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    add_subdirectory(tests)
    add_test(NAME SpaMtrixTests COMMAND SpaMtrixTests)
endif()