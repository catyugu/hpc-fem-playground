file(GLOB MPFEM_SOURCES "*.cpp")

add_library(mpfem STATIC ${MPFEM_SOURCES})
target_include_directories(mpfem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(mpfem PUBLIC cxx_std_17)
