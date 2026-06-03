add_library(affine AffineIntegrators.cpp AffineStressStrain.cpp)
target_include_directories(affine PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(affine PUBLIC mfem)
