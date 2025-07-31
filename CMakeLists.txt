cmake_minimum_required(VERSION 4.0.3)
project(EIT_cooling VERSION 0.1)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include path
include_directories(${CMAKE_SOURCE_DIR}/external/eigen)
include_directories(include)

# Core library sources
set(SOURCE_FILES
    src/chance_of_jump.cpp
    src/hamiltonian.cpp
    src/stochastic_schrodinger.cpp
)

# Header files (optional: for IDEs and clarity)
set(HEADER_FILES
    include/constants.hpp
    include/chance_of_jump.hpp
    include/hamiltonian.hpp
    include/stochastic_schrodinger.hpp
)

# Create static library
add_library(quantum_core STATIC ${SOURCE_FILES} ${HEADER_FILES})

# Example executable: fano_test
add_executable(fano_test examples/fano_test.cpp)
target_link_libraries(fano_test PRIVATE quantum_core)

# Optional: enable warnings
if (MSVC)
    target_compile_options(quantum_core PRIVATE /W4 /permissive-)
else()
    target_compile_options(quantum_core PRIVATE -Wall -Wextra -pedantic)
endif()

# Enable testing if desired
enable_testing()
add_subdirectory(tests)
