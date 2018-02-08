add_library(ACTARutils SHARED Utils.cc )
target_link_libraries(ACTARutils ${ROOT_LIBRARIES})


install(FILES Utils.h DESTINATION ${CMAKE_BINARY_DIR}/include)
