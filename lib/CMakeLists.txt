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
)

target_include_directories(SpaMtrix PUBLIC include)
install (TARGETS SpaMtrix DESTINATION bin)
install (FILES ${SPAMTRIX_LIB_INCLUDES} DESTINATION include)

