
link_directories(${SPAMTRIX_BINARY_DIR}/spamtrix) # libray path
add_executable(eigenvalues eigenvalues.cpp)   # add current executable
target_link_libraries(eigenvalues spamtrix)     # link library to to current executable
