
link_directories(${SPAMTRIX_BINARY_DIR}/spamtrix) # libray path
add_executable(pcg_test pcg_test.cpp)   # add current executable
target_link_libraries(pcg_test spamtrix)     # link library to to current executable
