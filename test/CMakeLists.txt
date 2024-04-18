# add unit tests
package_add_test(mqt-misim-test mqt-misim test_pkg.cpp)
package_add_test(mqt-misim-example mqt-misim main.cpp)

add_executable(mqt-misim-collect data_collect.cpp)
target_link_libraries(mqt-misim-collect PRIVATE mqt-misim)
set_target_properties(mqt-misim-collect PROPERTIES FOLDER tests)
