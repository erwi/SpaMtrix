
link_directories(${SPAMTRIX_BINARY_DIR}/spamtrix) # libray path
add_executable(input_output input_output.cpp)   # add current executable
target_link_libraries(input_output SpaMtrix)     # link library to to current executable
