# tests/CMakeLists.txt

add_executable(test_test test_test.cpp)
add_executable(test_chance_of_jump test_chance_of_jump.cpp)
add_executable(test_overlap_from_k test_overlap_from_k.cpp)
target_link_libraries(test_test PRIVATE quantum_core)
target_link_libraries(test_chance_of_jump PRIVATE quantum_core)
target_link_libraries(test_overlap_from_k PRIVATE quantum_core)

# Optionally integrate with CTest
enable_testing()
add_test(NAME TestTest COMMAND test_test)
add_test(NAME TestChanceOfJump COMMAND test_chance_of_jump)
add_test(NAME TestOverlapFromK COMMAND test_overlap_from_k)
