# MFEM Plugins Library

# Header-only library for plugins
add_library(mfem_plugins INTERFACE)

target_include_directories(mfem_plugins INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/plugins
)

# Note: mfem_plugins doesn't directly link to mfem here
# It will be linked through HPC_FEM_PLAYGROUND_LIB in parent CMakeLists.txt

# Install headers
install(DIRECTORY plugins/
    DESTINATION include/mfem_plugins
    FILES_MATCHING PATTERN "*.hpp"
)
