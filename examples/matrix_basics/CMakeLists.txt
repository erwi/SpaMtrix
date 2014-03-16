
link_directories(${SPAMTRIX_BINARY_DIR}/spamtrix) # libray path
add_executable(matrix_basics matrix_basics.cpp)   # add current executable
target_link_libraries(matrix_basics spamtrix)     # link library to to current executable
