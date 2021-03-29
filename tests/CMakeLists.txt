add_executable(SpaMtrixTests test.cpp)
target_include_directories(SpaMtrixTests
        PRIVATE
            ../extern/catch)

target_link_libraries(SpaMtrixTests SpaMtrix)