
link_directories(${SPAMTRIX_BINARY_DIR}/spamtrix) # libray path
add_executable(test6 test6.cpp)   # add current executable
target_link_libraries(test6 spamtrix) # link library to to current executable
