add_executable(test_id_range_parser test_id_range_parser.cpp)
target_link_libraries(test_id_range_parser PRIVATE mpfem)

add_executable(test_case_xml_summary test_case_xml_summary.cpp)
target_link_libraries(test_case_xml_summary PRIVATE mpfem)

add_executable(test_case_xml_reader_full test_case_xml_reader_full.cpp)
target_link_libraries(test_case_xml_reader_full PRIVATE mpfem)

add_executable(test_material_xml_reader test_material_xml_reader.cpp)
target_link_libraries(test_material_xml_reader PRIVATE mpfem)

add_executable(test_mfem_coupled_solver test_mfem_coupled_solver.cpp)
target_link_libraries(test_mfem_coupled_solver PRIVATE mpfem)

add_test(
    NAME test_id_range_parser
    COMMAND test_id_range_parser
)

add_test(
    NAME test_case_xml_summary
    COMMAND test_case_xml_summary ${PROJECT_SOURCE_DIR}/cases/busbar/case.xml
)

add_test(
    NAME test_case_xml_reader_full
    COMMAND test_case_xml_reader_full ${PROJECT_SOURCE_DIR}/cases/busbar/case.xml
)

add_test(
    NAME test_material_xml_reader
    COMMAND test_material_xml_reader ${PROJECT_SOURCE_DIR}/cases/busbar/material.xml
)

add_test(
    NAME test_mfem_coupled_solver
    COMMAND test_mfem_coupled_solver
            ${PROJECT_SOURCE_DIR}/cases/busbar/case.xml
            ${PROJECT_SOURCE_DIR}/cases/busbar/material.xml
            ${PROJECT_SOURCE_DIR}
)
