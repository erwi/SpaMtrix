# Create the libSpaMtrix.a static library
add_library(SpaMtrix ${SPAMTRIX_LIB_SRC})

install (TARGETS SpaMtrix DESTINATION bin)
install (FILES ${SPAMTRIX_LIB_INCLUDES} DESTINATION include)

