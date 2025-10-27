# Custom Methods of computing to extend MFEM functionality
# MFEM Plugins Library

# Header-only library for plugins
add_library(mfem_plugins INTERFACE)

target_include_directories(mfem_plugins INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/plugins
)

# Install headers
install(DIRECTORY plugins/
    DESTINATION include/mfem_plugins
    FILES_MATCHING PATTERN "*.hpp"
)
