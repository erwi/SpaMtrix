
link_directories(${SPAMTRIX_BINARY_DIR}/spamtrix) # libray path
add_executable(test5 test5.cpp)   # add current executable
target_link_libraries(test5 SpaMtrix)     # link library to to current executable
