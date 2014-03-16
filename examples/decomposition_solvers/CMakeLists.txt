# libray path
link_directories(${SPAMTRIX_BINARY_DIR}/spamtrix)
# add current executable
add_executable(decomposition_solvers decomposition_solvers.cpp)
# link library to to current executable
target_link_libraries(decomposition_solvers spamtrix)
