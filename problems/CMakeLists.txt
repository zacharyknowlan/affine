add_executable(FA FungArtery.cpp)
target_link_libraries(FA PRIVATE mfem affine)

add_executable(FL FungLigament.cpp)
target_link_libraries(FL PRIVATE mfem affine)

add_executable(HEA HyperElasticArtery.cpp)
target_link_libraries(HEA PRIVATE mfem affine)

add_executable(HEL HyperElasticLigament.cpp)
target_link_libraries(HEL PRIVATE mfem affine)

add_executable(FUE FungUniaxialExtension.cpp)
target_link_libraries(FUE PRIVATE mfem affine)

add_executable(FSI FungSingleInclusion.cpp)
target_link_libraries(FSI PRIVATE mfem affine)
