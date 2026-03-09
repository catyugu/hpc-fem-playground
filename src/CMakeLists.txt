file(GLOB_RECURSE MPFEM_SOURCES "*.cpp")

add_library(mpfem STATIC ${MPFEM_SOURCES})
target_include_directories(
	mpfem
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
		${CMAKE_CURRENT_SOURCE_DIR}/model
		${CMAKE_CURRENT_SOURCE_DIR}/io
		${CMAKE_CURRENT_SOURCE_DIR}/physics
		${CMAKE_CURRENT_SOURCE_DIR}/assembly
		${CMAKE_CURRENT_SOURCE_DIR}/solver
		${CMAKE_CURRENT_SOURCE_DIR}/coupling
		${CMAKE_CURRENT_SOURCE_DIR}/output
)
target_compile_features(mpfem PUBLIC cxx_std_17)
target_link_libraries(mpfem PUBLIC mfem)
