#====================
# BUILD THE LIBRARY
#====================
add_library(SpaMtrix
        src/cholesky.cpp
        src/cholincpreconditioner.cpp
        src/densematrix.cpp
        src/diagpreconditioner.cpp
        src/fleximatrix.cpp
        src/ircmatrix.cpp
        src/iterativesolvers.cpp
        src/lu.cpp
        src/luincpreconditioner.cpp
        src/matrixmaker.cpp
        src/powermethod.cpp
        src/reader.cpp
        src/spamtrix_blas.cpp
        src/tdmatrix.cpp
        src/vector.cpp
        src/writer.cpp
        include/spamtrix_exception.hpp
)
# public headers that should be automatically added to the includes search path
# of any other project that includes this as a sub-project
target_include_directories(SpaMtrix PUBLIC include)

#===================
# INSTALLATION
#===================
# List of headers that are copied to the installation destination.
set(INSTALL_HEADERS
        include/spamtrix_cholesky.hpp
        include/spamtrix_cholincpreconditioner.hpp
        include/spamtrix_densematrix.hpp
        include/spamtrix_diagpreconditioner.hpp
        include/spamtrix_fleximatrix.hpp
        include/spamtrix_ircmatrix.hpp
        include/spamtrix_iterativesolvers.hpp
        include/spamtrix_lu.hpp
        include/spamtrix_luincpreconditioner.hpp
        include/spamtrix_matrixmaker.hpp
        include/spamtrix_powermethod.hpp
        include/spamtrix_reader.hpp
        include/spamtrix_blas.hpp
        include/spamtrix_tdmatrix.hpp
        include/spamtrix_vector.hpp
        include/spamtrix_writer.hpp)

# set the "public headers" property.
set_target_properties(SpaMtrix PROPERTIES PUBLIC_HEADER "${INSTALL_HEADERS}")

# copy the library binary to the "bin" subdirectory and the public headers to
# the "include" subdirectory.
install(TARGETS SpaMtrix
        LIBRARY DESTINATION bin
        PUBLIC_HEADER DESTINATION include)
