
link_directories(${SPAMTRIX_BINARY_DIR}/spamtrix) # libray path
add_executable(test4 test4.cpp)   # add current executable
target_link_libraries(test4 spamtrix)     # link library to to current executable
