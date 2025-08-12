# tests/CMakeLists.txt

add_executable(test_test test_test.cpp)
add_executable(test_chance_of_jump test_chance_of_jump.cpp)
add_executable(test_overlap_from_k test_overlap_from_k.cpp)
add_executable(test_magnus_2 test_magnus2.cpp)
add_executable(test_wigner_symbol test_wigner_symbol.cpp)
add_executable(test_caf_structure test_caf_structure.cpp)
target_link_libraries(test_test PRIVATE quantum_core)
target_link_libraries(test_chance_of_jump PRIVATE quantum_core)
target_link_libraries(test_overlap_from_k PRIVATE quantum_core)
target_link_libraries(test_magnus_2 PRIVATE quantum_core)
target_link_libraries(test_wigner_symbol PRIVATE quantum_core)
target_link_libraries(test_caf_structure PRIVATE quantum_core)

# Optionally integrate with CTest
enable_testing()
add_test(NAME TestTest COMMAND test_test)
add_test(NAME TestChanceOfJump COMMAND test_chance_of_jump)
add_test(NAME TestOverlapFromK COMMAND test_overlap_from_k)
add_test(NAME TestMagnus2 COMMAND test_magnus_2)
add_test(NAME TestWignerSymbol COMMAND test_wigner_symbol)
add_test(NAME TestCafStructure COMMAND test_caf_structure)
