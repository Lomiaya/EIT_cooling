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
include_directories(include/CaF)

# Core library sources
set(SOURCE_FILES
    src/chance_of_jump.cpp
    src/hamiltonian.cpp
    src/stochastic_schrodinger.cpp
    src/stochastic_schrodinger_spin_motion.cpp
    src/atomic_states.cpp
    src/simulation.cpp
    src/CaF/wigner_symbol.cpp
    src/CaF/hunds_case_b.cpp
    src/CaF/state_hamiltonian.cpp
    src/caf_states.cpp
)

# Header files (optional: for IDEs and clarity)
set(HEADER_FILES
    include/constants.hpp
    include/chance_of_jump.hpp
    include/hamiltonian.hpp
    include/stochastic_schrodinger.hpp
    include/stochastic_schrodinger_spin_motion.hpp
    include/define_params.hpp
    include/atomic_states.hpp
    include/simulation.hpp
    include/CaF/half_integer.hpp
    include/CaF/hunds_case_b.hpp
    include/CaF/wigner_symbol.hpp
    include/CaF/state_hamiltonian.hpp
    include/caf_states.hpp
)

# Create static library
add_library(quantum_core STATIC ${SOURCE_FILES} ${HEADER_FILES})

# Example executable: fano_test
add_executable(fano_test examples/fano_test.cpp)
target_link_libraries(fano_test PRIVATE quantum_core)

# Example executable: simulation
add_executable(simulation examples/simulation.cpp)
target_link_libraries(simulation PRIVATE quantum_core)

# Example executable: simulation of CaF
add_executable(simulation_caf examples/simulation_caf.cpp)
target_link_libraries(simulation_caf PRIVATE quantum_core)

if(APPLE)
  message("is APPLE!")
  find_package(OpenMP REQUIRED)
  target_link_libraries(quantum_core PRIVATE OpenMP::OpenMP_CXX)
endif()

if (MSVC AND CMAKE_BUILD_TYPE STREQUAL "Release")
    message("is MSVC!")
    target_compile_options(quantum_core PRIVATE "/DNDEBUG")
    target_compile_options(quantum_core PRIVATE "/RTC-")
    target_compile_options(quantum_core PRIVATE "/Ox")
endif()

if (UNIX)
    message("is UNIX!")
    target_compile_options(quantum_core PRIVATE -DNDEBUG)
    target_compile_options(quantum_core PRIVATE -O3)
    find_package(Threads REQUIRED)
    target_link_libraries(quantum_core PRIVATE Threads::Threads)
endif()

# Enable testing if desired
# enable_testing()
# add_subdirectory(tests)
