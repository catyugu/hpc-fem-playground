file(GLOB_RECURSE EXAMPLE_SOURCES "*.cpp")
foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_SOURCE} NAME_WE)
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_SOURCE})
    target_include_directories(${EXAMPLE_NAME} PRIVATE 
        ${CMAKE_SOURCE_DIR}/src
    )
    target_link_libraries(${EXAMPLE_NAME} PRIVATE 
        HPC_FEM_PLAYGROUND_LIB
    )
endforeach()
