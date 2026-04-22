add_executable(SpaMtrixTests
        test_main.cpp
        test.cpp
        solver_test.cpp
        ircmatrix_test.cpp
        vector_test.cpp)
target_include_directories(SpaMtrixTests
        PRIVATE
            ../extern/catch)

target_link_libraries(SpaMtrixTests SpaMtrix)

if(BUILD_TESTING)
        set(SPAMTRIX_TEST_SOURCES
                ${CMAKE_CURRENT_SOURCE_DIR}/test.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/solver_test.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/ircmatrix_test.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/vector_test.cpp)

        # Re-run CMake automatically when a test source file changes so the
        # discovered Catch test list stays in sync with the implementation.
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${SPAMTRIX_TEST_SOURCES})

        function(spamtrix_register_catch_tests test_source)
                file(STRINGS "${test_source}" catch_test_lines REGEX "^[ \t]*TEST_CASE\\(\"[^\"]+\"")

                foreach(catch_test_line IN LISTS catch_test_lines)
                        string(REGEX REPLACE ".*TEST_CASE\\(\"([^\"]+)\".*" "\\1" catch_test_name "${catch_test_line}")
                        string(MAKE_C_IDENTIFIER "${catch_test_name}" catch_test_id)
                        get_filename_component(test_source_name "${test_source}" NAME_WE)
                        string(MAKE_C_IDENTIFIER "${test_source_name}_${catch_test_id}" ctest_name)

                        add_test(NAME "${ctest_name}" COMMAND SpaMtrixTests "${catch_test_name}")
                        set_tests_properties("${ctest_name}" PROPERTIES LABELS "${test_source_name}")
                endforeach()
        endfunction()

        foreach(test_source IN LISTS SPAMTRIX_TEST_SOURCES)
                spamtrix_register_catch_tests("${test_source}")
        endforeach()
endif()