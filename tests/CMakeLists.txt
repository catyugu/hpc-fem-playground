add_executable(test_id_range_parser test_id_range_parser.cpp)
target_link_libraries(test_id_range_parser PRIVATE mpfem)

add_executable(test_case_xml_summary test_case_xml_summary.cpp)
target_link_libraries(test_case_xml_summary PRIVATE mpfem)

add_test(
    NAME test_id_range_parser
    COMMAND test_id_range_parser
)

add_test(
    NAME test_case_xml_summary
    COMMAND test_case_xml_summary ${PROJECT_SOURCE_DIR}/cases/busbar/case.xml
)
