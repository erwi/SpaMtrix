
link_directories(${SPAMTRIX_BINARY_DIR}/spamtrix) # libray path
add_executable(tridiag_solver tridiag_solver.cpp)   # add current executable
target_link_libraries(tridiag_solver spamtrix)     # link library to to current executable
