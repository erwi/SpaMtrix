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
        add_test(NAME Vector_construction_initializes_all_entries_to_zero COMMAND SpaMtrixTests "Vector construction initializes all entries to zero")
        add_test(NAME Vector_copy_construction_creates_an_independent_copy COMMAND SpaMtrixTests "Vector copy construction creates an independent copy")
        add_test(NAME Vector_assignment_creates_an_independent_deep_copy COMMAND SpaMtrixTests "Vector assignment creates an independent deep copy")
        add_test(NAME Vector_scalar_assignment_sets_every_element COMMAND SpaMtrixTests "Vector scalar assignment sets every element")
        add_test(NAME Vector_element_access_supports_read_and_write COMMAND SpaMtrixTests "Vector element access supports read and write")
        add_test(NAME Vector_compound_arithmetic_with_scalars_and_vectors_works COMMAND SpaMtrixTests "Vector compound arithmetic with scalars and vectors works")
        add_test(NAME Vector_non_mutating_arithmetic_returns_new_vectors COMMAND SpaMtrixTests "Vector non-mutating arithmetic returns new vectors")
        add_test(NAME Vector_setAllValuesTo_norm_and_normalise_behave_correctly COMMAND SpaMtrixTests "Vector setAllValuesTo, norm, and normalise behave correctly")
        add_test(NAME Vector_resize_preserves_existing_values_and_zero_fills_new_elements COMMAND SpaMtrixTests "Vector resize preserves existing values and zero-fills new elements")

        add_test(NAME MatrixMaker_creates_IRCMatrix COMMAND SpaMtrixTests "MatrixMaker creates IRCMatrix")
        add_test(NAME Expand_block_along_matrix_diagonal COMMAND SpaMtrixTests "Expand block along matrix diagonal")
        add_test(NAME Expand_block_rows_and_columns COMMAND SpaMtrixTests "Expand block rows and columns")

        add_test(NAME Add_sparse_matrix_to_other_identical_sparsity_patterns COMMAND SpaMtrixTests "Add sparse matrix to other - identical sparsity patterns")
        add_test(NAME Add_sparse_matrix_to_other_different_sparsity_patterns COMMAND SpaMtrixTests "Add sparse matrix to other - different sparsity patterns")
        add_test(NAME Access_matrix_values_by_pointer COMMAND SpaMtrixTests "Access matrix values by pointer")

        add_test(NAME Solve_symmetric_positive_definite_5_point_Poisson_finite_difference_system COMMAND SpaMtrixTests "Solve symmetric positive definite 5-point Poisson finite difference system")
endif()